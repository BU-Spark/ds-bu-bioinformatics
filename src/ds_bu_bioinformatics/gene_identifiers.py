"""Cached, exact human identifier mappings for Entrez/Ensembl source matrices."""
import re
import urllib.parse

import numpy as np
import pandas as pd

from .sources import MYGENE_URL


def annotate(cache, feature_ids, *, scope):
    if scope not in {'entrezgene', 'ensembl.gene', 'symbol'}:
        raise ValueError('Unsupported source identifier namespace')
    pattern = {'entrezgene': r'\d+', 'ensembl.gene': r'ENSG\d+', 'symbol': r'[A-Za-z][A-Za-z0-9_.-]*'}[scope]
    queries = [re.sub(r'\.\d+$', '', x) if scope == 'ensembl.gene'
               else re.sub(r'_at$', '', x) if scope == 'entrezgene' else x
               for x in feature_ids]
    valid = sorted({q for q in queries if re.fullmatch(pattern, q)})
    by_query = {}
    for offset in range(0, len(valid), 1000):
        batch = valid[offset:offset + 1000]
        body = urllib.parse.urlencode({'q': ','.join(batch), 'scopes': scope,
            'fields': 'symbol,entrezgene,taxid,ensembl.gene', 'species': 'human', 'size': 100}).encode()
        hits = cache.json(f'mygene/{scope}/batch_{offset // 1000:04d}.json', MYGENE_URL, body)
        if not isinstance(hits, list) or {str(h.get('query', '')) for h in hits} != set(batch):
            raise ValueError('Identifier mapping response omits or adds source queries')
        for hit in hits:
            by_query.setdefault(str(hit['query']), []).append(hit)
    rows, candidate_rows = [], []
    for feature, query in zip(feature_ids, queries):
        identities, invalid = set(), False
        hits = by_query.get(query, [])
        for hit in hits:
            if hit.get('notfound'):
                continue
            gene, symbol = str(hit.get('entrezgene', '')), hit.get('symbol', '')
            if scope == 'entrezgene':
                exact = gene == query
            elif scope == 'ensembl.gene':
                ens = hit.get('ensembl', [])
                ens = [ens] if isinstance(ens, dict) else ens
                exact = any(e.get('gene') == query for e in ens)
            else:
                exact = symbol == query
            good = (exact and hit.get('taxid') == 9606 and gene.isdigit()
                    and isinstance(symbol, str) and bool(re.fullmatch(r'[^;/,\s]+', symbol)))
            invalid |= not good or int(hit.get('_total', 0)) >= 100
            if good:
                identities.add((gene, symbol))
            candidate_rows.append({'feature_id': feature, 'query': query, 'candidate_gene_id': gene,
                                   'candidate_symbol': symbol, 'exact_human_identity': good})
        status = 'mapped' if len(identities) == 1 and not invalid and len(hits) < 100 else 'invalid_or_ambiguous_identifier'
        gene, symbol = next(iter(identities)) if status == 'mapped' else ('', '')
        rows.append({'feature_id': feature, 'source_gene_id': feature, 'mapping_query': query,
                     'gene_id': gene, 'gene_symbol': symbol, 'mapping_status': status})
    frame = pd.DataFrame(rows)
    mapped = frame.loc[frame.mapping_status.eq('mapped'), ['gene_id', 'gene_symbol']].drop_duplicates()
    bad_ids = set(mapped.loc[mapped.gene_id.duplicated(False), 'gene_id'])
    bad_symbols = set(mapped.loc[mapped.gene_symbol.duplicated(False), 'gene_symbol'])
    conflict = frame.gene_id.isin(bad_ids) | frame.gene_symbol.isin(bad_symbols)
    frame.loc[conflict, 'mapping_status'] = 'gene_symbol_identity_conflict'
    frame['annotation_eligible'] = frame.mapping_status.eq('mapped')
    candidates = pd.DataFrame(candidate_rows, columns=['feature_id', 'query', 'candidate_gene_id', 'candidate_symbol', 'exact_human_identity'])
    return frame, candidates
