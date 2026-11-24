"""Nakagawa GSE10645: exact two-panel subject join and cached RefSeq annotation."""
import csv
import re

import numpy as np
import pandas as pd

from .annotation import fetch_mappings, resolve_mappings
from .cambridge import read_matrix
from .geo_cohort_pipeline import add_outcomes, decode_bcr, number

SERIES = 'GSE10645'
PANELS = ('GPL5858', 'GPL5873')
PLATFORM = ';'.join(PANELS)
ACCESSIONS = [SERIES, *PANELS]
STUDY = 'Nakagawa / Mayo (nested case-control; paired DASL panels)'
SNAPSHOT = 'nakagawa-v1'
SCALE = 'Deposited cyclic-loess/fastlo normalized gene-target signals; no added transformation'
MATRIX_URLS = {p: f'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/{SERIES}/matrix/{SERIES}-{p}_series_matrix.txt.gz' for p in PANELS}
PLATFORM_URLS = {p: f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={p}&targ=self&form=text&view=full' for p in PANELS}


def outcome(record):
    # Years of follow-up are not silently substituted for a separately verified
    # time to PCa death. A linked clinical supplement may provide that duration.
    return decode_bcr(record.get('geo__prostate_cancer_specific_death'), None, {'Yes': 1, 'No': 0})


ENDPOINTS = {'pcsm': outcome}


def panel_annotations(path, platform, expression):
    with path.open(encoding='utf-8') as handle:
        correct = False
        for line in handle:
            correct |= line.rstrip() == f'^PLATFORM = {platform}'
            if line.rstrip() == '!platform_table_begin':
                break
        else:
            raise ValueError(f'{platform}: platform table missing')
        reader = csv.reader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
        header = next(reader, [])
        if not correct or not {'ID', 'SequenceSource', 'GB_ACC'}.issubset(header) or len(set(header)) != len(header):
            raise ValueError(f'{platform}: unexpected DASL annotation schema')
        rows = []
        for row in reader:
            if row == ['!platform_table_end']:
                break
            if len(row) != len(header):
                raise ValueError('Malformed DASL annotation row')
            rows.append([row[header.index(c)].strip() for c in ('ID', 'SequenceSource', 'GB_ACC')])
        else:
            raise ValueError('Truncated DASL annotation')
    result = pd.DataFrame(rows, columns=['source_feature_id', 'sequence_source', 'refseq_original'])
    if result.empty or result.source_feature_id.eq('').any() or not result.source_feature_id.is_unique:
        raise ValueError('Missing/duplicate DASL source targets')
    if set(result.source_feature_id) != set(expression.index):
        raise ValueError('DASL expression/platform target sets differ; no silent intersection')
    result = result.set_index('source_feature_id').loc[expression.index.tolist()].reset_index()
    result['feature_id'] = platform + ':' + result.source_feature_id
    result['source_platform'] = platform
    result['refseq_query'] = result.refseq_original.str.replace(r'\.\d+$', '', regex=True)
    result['valid_refseq'] = (result.sequence_source.eq('RefSeq')
                              & result.refseq_query.str.fullmatch(r'(?:NM|NR|XM|XR)_\d+'))
    return result


def join_panels(panels):
    if set(panels) != set(PANELS):
        raise ValueError('Nakagawa requires both deposited panels')
    indexed, assay_rows, blocks = {}, [], []
    used_assays = set()
    for platform in PANELS:
        metadata, expression = panels[platform]
        meta = metadata.copy()
        if (meta.empty or not meta.platform_id.eq(platform).all()
                or meta.sample_id.tolist() != expression.columns.tolist()):
            raise ValueError('DASL platform/assay identity or order mismatch')
        subject = meta.title.str.extract(r'^Patient (\d+) Core Set(?: \([^()]+ panel\))?$', expand=False)
        if subject.isna().any() or not subject.is_unique or not meta.sample_id.is_unique:
            raise ValueError('DASL missing/duplicate patient label; no replicate averaging')
        if used_assays & set(meta.sample_id):
            raise ValueError('One GEO assay occurs on both panels')
        used_assays.update(meta.sample_id)
        meta['patient_id'] = subject
        meta['joint_sample_id'] = SERIES + ':patient:' + subject
        indexed[platform] = meta.set_index('patient_id', drop=False)
        assay_rows.append(meta.assign(assay_included=True))
    patient_ids = sorted(indexed[PANELS[0]].index)
    if set(patient_ids) != set(indexed[PANELS[1]].index):
        raise ValueError('DASL panels cover different subjects; obtain/reconcile missing panel instead of dropping patients')
    first, second = [indexed[p].loc[patient_ids].reset_index(drop=True) for p in PANELS]
    source_fields = sorted({c for frame in (first, second) for c in frame if c.startswith('geo__')})
    required = {'geo__prostate_cancer_specific_death', 'geo__years_of_followup_after_rrp', 'geo__case_control_group'}
    if not required.issubset(first) or not required.issubset(second):
        raise ValueError('DASL clinical status, group or follow-up source fields missing')
    for column in source_fields:
        if column not in first or column not in second:
            raise ValueError(f'Clinical source schema differs between panels: {column}')
        a = first[column].fillna('').astype(str).str.strip()
        b = second[column].fillna('').astype(str).str.strip()
        if not a.equals(b):
            raise ValueError(f'Conflicting panel clinical records: {column}; adjudicate source')
    clinical = first.copy()
    clinical['sample_id'] = clinical.joint_sample_id
    clinical['geo_accession'] = first.sample_id + ';' + second.sample_id
    clinical['platform_id'] = PLATFORM
    clinical['title'] = 'Patient ' + clinical.patient_id + ' paired DASL panels'
    for platform, frame in zip(PANELS, (first, second)):
        clinical[platform.lower() + '_assay_id'] = frame.sample_id
        expression = panels[platform][1].loc[:, frame.sample_id].copy()
        expression.columns = clinical.sample_id
        expression.index = pd.Index(platform + ':' + expression.index.astype(str), name='feature_id')
        blocks.append(expression)
    clinical = add_outcomes(clinical, outcome, 'pcsm')
    clinical['pcsm_clinical_source'] = 'GEO GSE10645 PCa-specific death label; separate event/censoring duration not verified'
    for source_key, target in [('geo__years_of_followup_after_rrp', 'source_followup_months_from_rrp'),
        ('geo__time_yrs_from_rrp_to_first_rising_psa', 'source_first_psa_months_from_rrp'),
        ('geo__time_yrs_from_rrp_to_second_rising_psa', 'source_second_psa_months_from_rrp'),
        ('geo__time_yrs_from_psa_recurrence_to_systemic_disease', 'source_systemic_interval_months_from_psa')]:
        if source_key in clinical:
            clinical[target] = pd.array([None if (v := number(x)) is None else 12 * v for x in clinical[source_key]], dtype='Float64')
    if not clinical.geo__case_control_group.isin(['Systemic', 'SYS', 'PSA', 'NED']).all():
        raise ValueError('Unexpected Nakagawa case-control group')
    samples = clinical[['sample_id', 'geo_accession', 'patient_id', 'platform_id', 'title',
                        'gpl5858_assay_id', 'gpl5873_assay_id']].copy()
    # The selected unit is a composite of two assays, not a claim that either assay
    # was independently selected as a second patient. Keep the full assay ledger too.
    audit = samples.assign(selected=True, selection_reason='paired_panels_exact_patient_join')
    return clinical, samples, audit, pd.concat(blocks), pd.concat(assay_rows, ignore_index=True)


def prepare(cache):
    panels, platforms, processing = {}, [], []
    for platform in PANELS:
        metadata, expression, methods = read_matrix(cache.fetch(f'{SERIES}-{platform}_series_matrix.txt.gz', MATRIX_URLS[platform]),
                                                     series=SERIES, platform=platform)
        panels[platform] = metadata, expression
        platforms.append(panel_annotations(cache.fetch(platform + '.soft.txt', PLATFORM_URLS[platform]), platform, expression))
        processing.extend(methods)
    clinical, samples, audit, source, assay_audit = join_panels(panels)
    platform_table = pd.concat(platforms, ignore_index=True)
    candidates = resolve_mappings(platform_table, fetch_mappings(cache, platform_table), source.index)
    collapsed = candidates.drop_duplicates('feature_id').set_index('feature_id').loc[source.index].reset_index()
    annotations = platform_table.set_index('feature_id').loc[source.index].reset_index()
    annotations['gene_id'] = collapsed.candidate_gene_id
    annotations['gene_symbol'] = collapsed.candidate_symbol
    annotations['mapping_status'] = collapsed.mapping_status
    annotations['annotation_eligible'] = collapsed.selected
    details = {
        'sample_selection': 'Both panels required for exactly the same title-derived patient labels; one composite sample per subject; all assay records audited',
        'patient_identity': 'Exact source Patient NNN Core Set label; clinical fields must agree across both panels',
        'annotation': 'Cached human RefSeq-to-Entrez/symbol mapping with verified returned RefSeq; no guessing from GI/Affymetrix target labels',
        'endpoints': {'pcsm': {'definition': 'Deposited prostate-cancer-specific death class', 'unit': 'months',
            'duration_available': False, 'missing_policy': 'PCSM times remain missing until a linked clinical supplement verifies event/censoring dates'}},
        'study_source': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC2565588/',
        'limitations': ['Complete-source acquisition/counts not yet verified in this environment',
            'Deposited GPL5858/GPL5873 gene-target signals already average probes; original raw probe preservation is not claimed',
            'Non-RefSeq GenBank and Affymetrix target annotations are excluded and audited; broader mapping requires explicit validation',
            'PSA control group includes later systemic events: never equate group label with final systemic-event status',
            'Years of follow-up are retained and converted separately; not automatically substituted for PCSM duration',
            'Nested case-control sampling and targeted gene coverage limit general-population survival interpretation']}
    return clinical, samples, audit, source, annotations, {'assay_linkage.csv': assay_audit,
        'gene_mapping_candidates.csv.gz': candidates}, sorted(set(processing)), details


def summarize(clinical, audit, mapping, extra):
    return {'assays': len(extra['assay_linkage.csv']), 'case_control_groups': clinical.geo__case_control_group.value_counts().to_dict(),
            'survival_ready': bool(clinical.pcsm_eligible.any()),
            'source_feature_policy': 'Preserve separate deposited gene-target rows, namespaced by panel; no new aggregation'}
