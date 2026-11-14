"""Long GSE54460: lossless FPKM rows, explicit assay linkage and patient selection."""
import csv
import re

import numpy as np
import pandas as pd

from .cambridge import read_matrix
from .expression import open_text
from .geo_cohort_pipeline import MISSING, add_outcomes, decode_bcr, number

SERIES = 'GSE54460'
PLATFORM = 'GPL11154'
STUDY = 'Long / FFPE RNA-seq'
SNAPSHOT = 'long-v1'
SCALE = 'Deposited untransformed FPKM; not counts, TPM or log2 expression'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54460/matrix/GSE54460_series_matrix.txt.gz'
FPKM_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54460/suppl/GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt.gz'


def outcome(record):
    status = record.get('geo__bcr')
    event_time = number(record.get('geo__months_to_bcr'))
    followup = number(record.get('geo__months_total_f_u'))
    duration = event_time if status == '1' else followup if status == '0' else None
    event, months, reason = decode_bcr(status, duration, {'1': 1, '0': 0})
    if event == 1 and months is not None and followup is not None and months > followup + 1e-8:
        reason = 'event_after_recorded_followup_requires_review'
    return event, months, reason


def read_fpkm(path):
    with open_text(path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, [])
        if len(header) < 3 or header[:2] != ['', '']:
            raise ValueError('Unexpected Long FPKM sample header')
        ids = header[2:]
        if any(not sid for sid in ids) or len(set(ids)) != len(ids):
            raise ValueError('Duplicate/empty FPKM assay identifiers')
        fields, values, annotations, bad = {}, [], [], []
        table = False
        for line_number, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError('Malformed Long FPKM row')
            if not table:
                if row[:2] == ['Entrez', 'SYMBOL']:
                    table = True
                elif any(row):
                    key = 'fpkm__' + re.sub('[^a-z0-9]+', '_', row[0].strip().lower()).strip('_')
                    if not row[0] or row[1] or key in fields:
                        raise ValueError('Invalid/duplicate embedded FPKM clinical header')
                    fields[key] = row[2:]
                continue
            if not row[0].isdigit() or not row[1]:
                raise ValueError('Missing/malformed deposited FPKM gene annotation')
            feature = f'{SERIES}:line:{line_number}'
            annotations.append({'feature_id': feature, 'gene_id': row[0], 'gene_symbol': row[1], 'source_row': line_number})
            parsed = []
            for assay, token in zip(ids, row[2:]):
                if token.strip().upper() in MISSING | {'#NUM!'}:
                    parsed.append(np.nan)
                    bad.append({'feature_id': feature, 'expression_sample_id': assay, 'source_token': token,
                                'reason': 'spreadsheet_numeric_error' if token == '#NUM!' else 'missing_source_value'})
                else:
                    value = float(token)
                    if not np.isfinite(value) or value < 0:
                        raise ValueError('Invalid/nonfinite or negative FPKM value')
                    parsed.append(value)
            values.append(parsed)
    if not table or not annotations or 'fpkm__external_id' not in fields:
        raise ValueError('Incomplete Long FPKM table')
    genes = pd.DataFrame(annotations)
    genes['gene_id_repeated_in_source'] = genes.gene_id.duplicated(False)
    genes['mapping_status'] = 'deposited_annotation_not_revalidated'
    genes['annotation_eligible'] = True
    source = pd.DataFrame(values, index=pd.Index(genes.feature_id, name='feature_id'), columns=ids, dtype=float)
    assay_metadata = pd.DataFrame({'expression_sample_id': ids, **fields})
    if assay_metadata.fpkm__external_id.isna().any() or assay_metadata.fpkm__external_id.eq('').any():
        raise ValueError('FPKM external patient IDs absent')
    errors = pd.DataFrame(bad, columns=['feature_id', 'expression_sample_id', 'source_token', 'reason'])
    return source, genes, assay_metadata, errors


def link_and_select(metadata, assay_metadata):
    required = {'geo__cgsr_id', 'geo__external_id', 'geo__no_of_uniquely_mapped_reads', 'geo__tissue', 'geo__bcr'}
    if not required.issubset(metadata) or not metadata.geo__cgsr_id.is_unique or not metadata.sample_id.is_unique:
        raise ValueError('Missing/duplicate Long GEO assay identity')
    if not metadata.geo__tissue.eq('prostate').all() or not metadata.title.str.startswith('Prostate_').all():
        raise ValueError('Unexpected Long tissue/title; primary-cohort selection needs review')
    for r in metadata.to_dict('records'):
        if r['title'] not in {'Prostate_' + r['geo__external_id'], 'Prostate_' + r['geo__external_id'] + '-2'}:
            raise ValueError('GEO title disagrees with external patient ID')
    if len(metadata) != len(assay_metadata) or not assay_metadata.expression_sample_id.is_unique:
        raise ValueError('Long expression and GEO assay counts/identities differ')
    by_cgsr = metadata.set_index('geo__cgsr_id')
    matches, methods, used = {}, {}, set()
    # Resolve exact assay keys first, then uniquely remaining external-ID aliases.
    for r in assay_metadata.to_dict('records'):
        sid = r['expression_sample_id']
        if sid in by_cgsr.index:
            row = by_cgsr.loc[sid]
            if row.geo__external_id != r['fpkm__external_id']:
                raise ValueError('Exact assay identity disagrees with external patient')
            matches[sid] = row.sample_id
            methods[sid] = 'exact_cgsr_id'
            used.add(row.sample_id)
    for r in assay_metadata.to_dict('records'):
        sid = r['expression_sample_id']
        if sid in matches:
            continue
        candidates = metadata.loc[metadata.geo__external_id.eq(r['fpkm__external_id']) & ~metadata.sample_id.isin(used)]
        if len(candidates) != 1:
            raise ValueError('Unresolved or ambiguous FPKM external-ID alias: ' + sid)
        matches[sid] = candidates.sample_id.iloc[0]
        methods[sid] = 'unique_remaining_external_id'
        used.add(matches[sid])
    if len(used) != len(metadata):
        raise ValueError('Long assay linkage is not a bijection')
    linked = assay_metadata.copy()
    linked['sample_id'] = linked.expression_sample_id.map(matches)
    linked['expression_link_method'] = linked.expression_sample_id.map(methods)
    audit = metadata.merge(linked, on='sample_id', how='left', validate='one_to_one', sort=False)
    audit['patient_id'] = audit.geo__external_id
    # Assay quality is an outcome-independent selection criterion. No means.
    quality = audit.geo__no_of_uniquely_mapped_reads.str.replace(',', '', regex=False)
    if not quality.str.fullmatch(r'\d+').all():
        raise ValueError('Malformed/missing uniquely mapped read count')
    audit['uniquely_mapped_reads'] = quality.astype('int64')
    if not audit.uniquely_mapped_reads.gt(0).all():
        raise ValueError('Nonpositive uniquely mapped read count')
    for patient, group in audit.groupby('patient_id', sort=False):
        statuses = {x for x in group.geo__bcr if str(x).strip().upper() not in MISSING}
        if not statuses.issubset({'0', '1'}) or len(statuses) > 1:
            raise ValueError('Conflicting or unknown BCR statuses for patient ' + patient)
    ranked = audit.sort_values(['patient_id', 'uniquely_mapped_reads', 'sample_id'], ascending=[True, False, True], kind='stable')
    winners = set(ranked.drop_duplicates('patient_id').sample_id)
    audit['selected'] = audit.sample_id.isin(winners)
    audit['selection_reason'] = np.where(audit.selected, 'primary_selected_highest_mapped_reads', 'replicate_lower_priority')
    clinical = add_outcomes(audit.loc[audit.selected].reset_index(drop=True), outcome)
    clinical['bcr_clinical_source'] = 'GEO GSE54460 sample characteristics; no mixing with FPKM embedded clinical times'
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title', 'expression_sample_id', 'expression_link_method', 'uniquely_mapped_reads']].copy()
    return clinical, samples, audit


def comparison_table(audit):
    result = audit[['sample_id', 'patient_id', 'expression_sample_id', 'expression_link_method', 'selected']].copy()
    for suffix, geo, fpkm in [('status', 'geo__bcr', 'fpkm__bcr'),
                              ('event_months', 'geo__months_to_bcr', 'fpkm__months_to_bcr'),
                              ('followup_months', 'geo__months_total_f_u', 'fpkm__months_total_f_u')]:
        if geo not in audit or fpkm not in audit:
            raise ValueError('Clinical comparison field absent: ' + suffix)
        result['geo_' + suffix] = [number(x) for x in audit[geo]]
        result['fpkm_' + suffix] = [number(x) for x in audit[fpkm]]
        result[suffix + '_agrees'] = np.isclose(result['geo_' + suffix].to_numpy(dtype=float), result['fpkm_' + suffix].to_numpy(dtype=float), atol=1e-8, rtol=0, equal_nan=True)
    return result


def prepare(cache):
    meta, empty, processing = read_matrix(cache.fetch('GSE54460_series_matrix.txt.gz', MATRIX_URL), series=SERIES, platform=PLATFORM, allow_empty=True)
    if not empty.empty:
        raise ValueError('Unexpected expression in Long metadata-only Series Matrix')
    raw, annotations, assay_meta, errors = read_fpkm(cache.fetch('GSE54460_FPKM.txt.gz', FPKM_URL))
    clinical, samples, audit = link_and_select(meta, assay_meta)
    order = samples.expression_sample_id.tolist()
    source = raw.loc[:, order].copy()
    source.columns = samples.sample_id.tolist()
    source.index.name = 'feature_id'
    errors['sample_id'] = errors.expression_sample_id.map(audit.set_index('expression_sample_id').sample_id)
    errors['sample_selected'] = errors.sample_id.isin(samples.sample_id)
    comparisons = comparison_table(audit)
    details = {'sample_selection': 'Source FFPE radical-prostatectomy cohort; GEO prostate tissue and Prostate_* title required. One assay/patient: greatest GEO uniquely mapped reads, then lowest GSM tie-break; independent of outcomes',
        'patient_identity': 'Exact matching GEO cgsr ID plus external patient ID where possible; otherwise unique remaining GEO external-ID record, requiring a full assay bijection',
        'clinical_source': 'GEO status, event time and follow-up consistently for every patient; embedded FPKM clinical fields retained for comparison, never substituted',
        'annotation': 'Original FPKM Entrez and symbol labels preserved without current annotation validation. Repeated Entrez IDs and differing labels are retained as separate source rows',
        'feature_key': 'GSE54460:line:<original 1-based physical source line>; Entrez is an annotation, not a unique row key',
        'invalid_expression': '#NUM! and recognized missing tokens audited; exclude entire feature if affected in any selected assay; no zero fill or averaging',
        'zero_event_times': 'Source 0-month events retained but ineligible; never add pseudotime or infer missing dates',
        'study_source': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4058362/',
        'limitations': ['GEO and embedded FPKM clinical times differ; authoritative adjudication remains open',
            'Deposited gene labels/Entrez associations may be outdated or inconsistent',
            'No raw read reprocessing or expression outlier exclusion; deposited FPKM retained unchanged']}
    extra = {'clinical_source_comparison.csv': comparisons, 'expression_value_audit.csv': errors}
    return clinical, samples, audit, source, annotations, extra, processing, details


def summarize(clinical, audit, mapping, extra):
    comparison = extra['clinical_source_comparison.csv']
    return {'assay_link_methods': audit.expression_link_method.value_counts().to_dict(),
        'replicate_assays_excluded': int((~audit.selected).sum()),
        'patients_with_multiple_assays': int(audit.groupby('patient_id').size().gt(1).sum()),
        'clinical_source_disagreements_all_assays': {x: int((~comparison[x + '_agrees']).sum()) for x in ['status', 'event_months', 'followup_months']},
        'clinical_source_disagreements_selected_assays': {x: int((~comparison.loc[comparison.selected, x + '_agrees']).sum()) for x in ['status', 'event_months', 'followup_months']},
        'expression_invalid_source_cells': len(extra['expression_value_audit.csv']),
        'source_rows_sharing_entrez_id': int(mapping.gene_id_repeated_in_source.sum())}
