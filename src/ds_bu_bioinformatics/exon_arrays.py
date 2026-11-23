"""Cohort-specific GEO acquisition for the four Decipher exon-array studies."""
import csv
import re

import pandas as pd

from .annotation import fetch_mappings, resolve_mappings
from .geo_cohort_pipeline import add_outcomes, decode_bcr
from .geo_io import align_columns, read_geo, read_tsv


def source_annotations(cache, platform, expression):
    url = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={platform}&targ=self&form=text&view=full'
    path = cache.fetch(platform + '.soft.txt', url)
    rows, seen = [], set()
    wanted = set(expression.index)
    with path.open(encoding='utf-8') as handle:
        correct = False
        for line in handle:
            correct |= line.rstrip() == f'^PLATFORM = {platform}'
            if line.rstrip() == '!platform_table_begin':
                break
        else:
            raise ValueError('Affymetrix platform table missing')
        reader = csv.reader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
        header = next(reader, [])
        if not correct or not {'ID', 'gene_assignment'}.issubset(header) or len(set(header)) != len(header):
            raise ValueError('Unexpected Affymetrix gene_assignment schema')
        id_position, assignment_position = header.index('ID'), header.index('gene_assignment')
        for row in reader:
            if row == ['!platform_table_end']:
                break
            if len(row) != len(header) or not row[id_position] or row[id_position] in seen:
                raise ValueError('Malformed/duplicate Affymetrix annotation')
            feature = row[id_position]
            seen.add(feature)
            if feature not in wanted:
                continue
            assignment = row[assignment_position]
            symbols, refseqs, stable_ids, malformed = set(), set(), set(), False
            if assignment not in {'', '---'}:
                for hit in assignment.split(' /// '):
                    fields = [v.strip() for v in hit.split(' // ')]
                    if len(fields) < 2:
                        malformed = True
                        continue
                    symbols.add(fields[1])
                    if len(fields) >= 5 and fields[4].isdigit():
                        stable_ids.add(fields[4])
                    query = re.sub(r'\.\d+$', '', fields[0])
                    if re.fullmatch(r'(?:NM|NR|XM|XR)_\d+', query):
                        refseqs.add(query)
            valid = len(symbols) == 1 and len(stable_ids) <= 1 and bool(refseqs) and not malformed
            rows.append({'feature_id': feature, 'refseq_original': ';'.join(sorted(refseqs)),
                         'refseq_query': sorted(refseqs) if valid else [''],
                         'valid_refseq': valid, 'source_gene_assignment': assignment,
                         'source_gene_id': next(iter(stable_ids)) if len(stable_ids) == 1 else '',
                         'source_gene_symbol': next(iter(symbols)) if len(symbols) == 1 else ''})
        else:
            raise ValueError('Truncated Affymetrix platform table')
    if not wanted.issubset(seen):
        raise ValueError('Expression feature absent from original platform; no silent intersection')
    table = pd.DataFrame(rows).set_index('feature_id').loc[expression.index].reset_index()
    expanded = table.explode('refseq_query').reset_index(drop=True)
    candidates = resolve_mappings(expanded, fetch_mappings(cache, expanded), expression.index)
    grouped = {feature: frame for feature, frame in candidates.groupby('feature_id', sort=False)}
    outcomes = []
    for row in table.itertuples():
        hits = grouped[row.feature_id]
        selected = hits.loc[hits.selected]
        identities = set(zip(selected.candidate_gene_id, selected.candidate_symbol))
        conflicting = (hits.mapping_status.isin(['ambiguous', 'invalid_candidate', 'refseq_mismatch',
                                                'gene_symbol_identity_conflict']).any()
                       or len(identities) > 1)
        gene, symbol = next(iter(identities)) if len(identities) == 1 else ('', '')
        conflicting |= bool(symbol and symbol != row.source_gene_symbol)
        conflicting |= bool(gene and row.source_gene_id and gene != row.source_gene_id)
        status = ('source_mapping_identity_conflict' if conflicting else 'mapped' if gene
                  else 'unmapped' if row.valid_refseq else 'invalid_source_assignment')
        outcomes.append((gene if status == 'mapped' else '', symbol if status == 'mapped' else '', status, status == 'mapped'))
    table['refseq_query'] = table.refseq_query.map(';'.join)
    table[['gene_id', 'gene_symbol', 'mapping_status', 'annotation_eligible']] = pd.DataFrame(outcomes, index=table.index)
    decisions = table.set_index('feature_id').mapping_status
    candidates['feature_mapping_status'] = candidates.feature_id.map(decisions)
    candidates['selected'] &= candidates.feature_mapping_status.eq('mapped')
    return table, candidates


class ExonCohort:
    PLATFORM = 'GPL5188'
    SCALE = 'Deposited normalized PSR-level expression; no added transformation or probe aggregation'

    def __init__(self, key, series, study, method, *, supplementary=None, postoperative_rt=False):
        self.key, self.SERIES, self.STUDY = key, series, study
        self.SNAPSHOT = key + '-v1'
        bucket = series[:-3] + 'nnn'
        base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{series}'
        self.MATRIX_URL = base + f'/matrix/{series}_series_matrix.txt.gz'
        self.SUPPLEMENT_URL = base + '/suppl/' + supplementary if supplementary else None
        self.supplementary = supplementary
        self.method, self.postoperative_rt = method, postoperative_rt
        self.ENDPOINTS = {'metastasis': self.outcome}
        if postoperative_rt:
            self.ENDPOINTS['bcr'] = self.bcr_outcome

    def outcome(self, record):
        if self.SERIES != 'GSE46691':
            return None, None, 'missing_status'
        status = record.get('geo__metastatic_event')
        return decode_bcr(status, None, {'0': 0, '1': 1})

    def bcr_outcome(self, record):
        # No public per-sample BCR field verified for this series. Do not recode
        # Gleason, PSA or treatment labels into recurrence outcomes.
        return None, None, 'missing_status'

    def clinical_tables(self, metadata):
        audit = metadata.copy()
        if (audit.empty or not audit.platform_id.eq(self.PLATFORM).all()
                or not audit.title.str.fullmatch(r'prostate_cancer_primary_\d+(?: \[[^\]]+\])?').all()
                or not audit.sample_id.is_unique):
            raise ValueError(f'{self.STUDY}: unexpected primary specimen identity/title')
        # Cohort-qualified source study subject label, never merged with another Mayo cohort.
        audit['patient_id'] = self.SERIES + ':' + audit.title.str.extract(r'^prostate_cancer_primary_(\d+)', expand=False)
        if not audit.patient_id.is_unique:
            raise ValueError('Duplicate study subject; no arbitrary replicate merging')
        audit['selected'] = True
        audit['selection_reason'] = 'primary_prostatectomy_tumor'
        clinical = audit.reset_index(drop=True)
        for prefix, decoder in self.ENDPOINTS.items():
            clinical = add_outcomes(clinical, decoder, prefix)
            clinical[f'{prefix}_clinical_source'] = self.SERIES + ' GEO source labels; individual endpoint durations unverified'
        if self.postoperative_rt:
            clinical['treatment_context'] = 'Study cohort received postoperative radiation; not per-patient treatment adjudication'
        samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title']].copy()
        return clinical, samples, audit

    def prepare(self, cache):
        meta, expression, processing = read_geo(cache.fetch(self.SERIES + '_series_matrix.txt.gz', self.MATRIX_URL),
                                                series=self.SERIES, platforms={self.PLATFORM})
        extra = {}
        if expression.empty:
            if not self.SUPPLEMENT_URL:
                raise ValueError(f'{self.SERIES}: source matrix empty; inspect the deposited processed files')
            expression = read_tsv(cache.fetch(self.supplementary, self.SUPPLEMENT_URL))
            expression, linkage = align_columns(expression, meta,
                lambda r: [r['sample_id'], r.get('source_description', '')])
            extra['expression_linkage.csv'] = linkage
        if set(expression.columns) != set(meta.sample_id):
            raise ValueError('Expression and GEO source assay sets differ')
        clinical, samples, audit = self.clinical_tables(meta)
        source = expression.loc[:, samples.sample_id]
        annotations, candidates = source_annotations(cache, self.PLATFORM, source)
        extra['gene_mapping_candidates.csv.gz'] = candidates
        details = {
            'sample_selection': 'Require primary-tumor title and expected platform; retain all unique study-qualified subjects independently of outcomes',
            'patient_identity': 'Source title-derived within-study label; overlap across Mayo/Decipher studies is not resolved',
            'annotation': 'Original PSR IDs retained; single-symbol source assignments anchored by cached, exact human RefSeq mapping; no aggregation',
            'expression_kind': 'log_microarray',
            'endpoints': {p: {'definition': 'Source metastatic event' if p == 'metastasis' else 'Biochemical recurrence',
                'duration_available': False, 'unit': 'months', 'missing_policy': 'Missing status/time remains missing; verified clinical supplement can add individual durations'} for p in self.ENDPOINTS},
            'study_source': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={self.SERIES}',
            'limitations': ['Full-source run and counts remain unverified; synthetic integration tests only',
                'Missing clinical follow-up does not block expression/PCA; survival fitting requires verified individual event/censoring times',
                'Source processed values (' + self.method + ') retained; no raw CEL QC/reprocessing',
                'Large PSR matrices require substantial memory/disk; numeric conversion is batched',
                'Case-cohort sampling, treatment context and patient overlap need analysis-specific review']}
        return clinical, samples, audit, source, annotations, extra, processing, details

    def summarize(self, clinical, audit, mapping, extra):
        return {'survival_ready': any(clinical[f'{p}_eligible'].any() for p in self.ENDPOINTS)}


ERHO = ExonCohort('erho', 'GSE46691', 'Erho / Mayo', 'fRMA', supplementary='GSE46691_quantile_normalized.txt.gz')
KARNES = ExonCohort('karnes', 'GSE62116', 'Karnes / Mayo', 'fRMA')
DEN = ExonCohort('den', 'GSE72291', 'Den / Thomas Jefferson', 'deposited normalized signal', postoperative_rt=True)
KLEIN = ExonCohort('klein', 'GSE62667', 'Klein / Cleveland', 'SCAN')
