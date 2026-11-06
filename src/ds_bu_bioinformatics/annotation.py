"""Preserve complete MyGene candidates and separate unambiguously mapped features."""
from __future__ import annotations
import csv
import re
import urllib.parse
import pandas as pd
from .sources import MYGENE_URL


def parse_platform(path):
    rows = []
    with path.open(encoding='utf-8') as handle:
        in_table = False
        ended = False
        header = None
        platform = False
        for line in handle:
            line = line.rstrip('\r\n')
            if line == '^PLATFORM = GPL10264':
                platform = True
            if line == '!platform_table_begin':
                in_table = True
                continue
            if line == '!platform_table_end':
                ended = True
                break
            if not in_table:
                continue
            parts = next(csv.reader([line], delimiter='\t'))
            if header is None:
                header = parts
                if not {'ID', 'GB_ACC'}.issubset(header):
                    raise ValueError('Unexpected GPL10264 annotation columns')
                continue
            if len(parts) != len(header):
                raise ValueError('Malformed GPL annotation row')
            rows.append((parts[header.index('ID')], parts[header.index('GB_ACC')]))
    result = pd.DataFrame(rows, columns=['feature_id', 'refseq_original'])
    if not platform or not ended or result.empty or result.feature_id.duplicated().any() or result.feature_id.eq('').any():
        raise ValueError('Wrong/truncated/duplicate GPL10264 annotation')
    result['refseq_query'] = result.refseq_original.str.replace(r'\.\d+$', '', regex=True)
    result['valid_refseq'] = result.refseq_query.str.fullmatch(r'(?:NM|NR|XM|XR)_\d+')
    return result


def fetch_mappings(cache, platform, batch_size=1000):
    queries = sorted(platform.loc[platform.valid_refseq, 'refseq_query'].unique())
    payload = []
    for start in range(0, len(queries), batch_size):
        batch = queries[start:start + batch_size]
        body = urllib.parse.urlencode({'q': ','.join(batch), 'scopes': 'refseq',
            'fields': 'symbol,entrezgene,taxid,refseq', 'species': 'human', 'size': 100}).encode()
        records = cache.json(f'mygene/batch_{start // batch_size:04d}.json', MYGENE_URL, body)
        if not isinstance(records, list) or any(str(r.get('query', '')) not in batch for r in records):
            raise ValueError('Unexpected MyGene query response')
        if set(str(r.get('query')) for r in records) != set(batch):
            raise ValueError('MyGene omitted query records')
        payload.extend(records)
    return payload


def resolve_mappings(platform, payload, expression_ids):
    if set(expression_ids) != set(platform.feature_id):
        raise ValueError('Expression/GPL feature sets differ; no silent inner join')
    by_query = {}
    for hit in payload:
        by_query.setdefault(str(hit['query']), []).append(hit)
    rows = []
    for feature in platform.itertuples(index=False):
        hits = by_query.get(feature.refseq_query, []) if feature.valid_refseq else []
        candidates = [hit for hit in hits if not hit.get('notfound')]
        identities = set()
        valid = True
        truncated = len(candidates) >= 100
        refseq_mismatch = False
        for hit in candidates:
            truncated |= int(hit.get('_total', 0)) >= 100
            gene_id, symbol = hit.get('entrezgene'), hit.get('symbol')
            rna = hit.get('refseq', {}).get('rna', [])
            rna = [rna] if isinstance(rna, str) else rna
            refseq_mismatch |= feature.refseq_query not in {re.sub(r'\.\d+$', '', x) for x in rna}
            if hit.get('taxid') != 9606 or not str(gene_id).isdigit() or not isinstance(symbol, str) or not symbol:
                valid = False
            else:
                identities.add((str(gene_id), symbol))
        status = ('invalid_refseq' if not feature.valid_refseq else 'unmapped' if not candidates else
                  'ambiguous' if len(identities) > 1 or truncated else 'invalid_candidate' if not valid or not identities else 'refseq_mismatch' if refseq_mismatch else 'mapped')
        for hit in candidates or [{}]:
            rows.append({'feature_id': feature.feature_id, 'refseq_original': feature.refseq_original,
                         'refseq_query': feature.refseq_query, 'candidate_gene_id': str(hit.get('entrezgene', '')),
                         'candidate_symbol': hit.get('symbol', ''), 'candidate_taxid': hit.get('taxid'),
                         'candidate_score': hit.get('_score'),
                         'candidate_refseq_rna': ';'.join(sorted([hit.get('refseq', {}).get('rna')] if isinstance(hit.get('refseq', {}).get('rna'), str) else hit.get('refseq', {}).get('rna', []))),
                         'mapping_status': status,
                         'selected': status == 'mapped', 'reason': status,
                         'mapping_source': 'MyGene.info v3 cached response'})
    mapping = pd.DataFrame(rows)
    # The same symbol must not label distinct stable gene identities; conflicting symbols
    # for one identity are likewise quarantined rather than merged by text label.
    selected = mapping.loc[mapping.selected, ['candidate_gene_id', 'candidate_symbol']].drop_duplicates()
    bad_symbols = set(selected.loc[selected.candidate_symbol.duplicated(False), 'candidate_symbol'])
    bad_ids = set(selected.loc[selected.candidate_gene_id.duplicated(False), 'candidate_gene_id'])
    conflict = mapping.selected & (mapping.candidate_symbol.isin(bad_symbols) | mapping.candidate_gene_id.isin(bad_ids))
    mapping.loc[conflict, ['mapping_status', 'reason']] = 'gene_symbol_identity_conflict'
    mapping.loc[conflict, 'selected'] = False
    return mapping.sort_values(['feature_id', 'candidate_gene_id', 'candidate_symbol'], kind='stable').reset_index(drop=True)


def select_features(expression, mapping):
    chosen = mapping[mapping.selected].drop_duplicates(['feature_id', 'candidate_gene_id', 'candidate_symbol'])
    if chosen.empty or chosen.feature_id.duplicated().any():
        raise ValueError('No mapped features or multiple selected mappings per feature')
    # Preserve source order and feature identity, even when gene annotations repeat.
    chosen = chosen.set_index('feature_id').loc[expression.index[expression.index.isin(chosen.feature_id)]]
    selected = expression.loc[chosen.index].copy()
    selected.index.name = 'feature_id'
    genes = chosen[['candidate_gene_id', 'candidate_symbol', 'refseq_original', 'refseq_query']].rename(
        columns={'candidate_gene_id': 'gene_id', 'candidate_symbol': 'gene_symbol'}).reset_index()
    return selected, genes
