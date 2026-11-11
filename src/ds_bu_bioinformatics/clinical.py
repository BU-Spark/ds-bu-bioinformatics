"""cBioPortal clinical snapshots and explicit GEO sample/patient matching."""
from __future__ import annotations
import json
import math
import urllib.parse
import pandas as pd
from .sources import CBIO_BASE

DFS_CODES = {'0:DiseaseFree': 0, '1:Recurred': 1, '1:Recurred/Progressed': 1}


def paged(cache, endpoint, prefix, params=None, page_size=1000):
    records = []
    fingerprints = set()
    expected = None
    for page in range(10000):
        query = urllib.parse.urlencode({**(params or {}), 'pageSize': page_size, 'pageNumber': page})
        name = f'{prefix}_{page:04d}.json'
        payload = cache.json(name, f'{CBIO_BASE}/{endpoint}?{query}')
        if not isinstance(payload, list):
            raise ValueError(f'Expected list for {endpoint}')
        entry = cache.manifest['sources'][name]
        total = entry['response_metadata'].get('total-count')
        if total is not None:
            if expected is not None and expected != int(total):
                raise ValueError('Clinical source changed during pagination')
            expected = int(total)
        fingerprint = json.dumps(payload, sort_keys=True)
        if payload and fingerprint in fingerprints:
            raise ValueError('Repeated API page; pagination is not progressing')
        fingerprints.add(fingerprint)
        records.extend(payload)
        if not payload:
            if expected is not None and len(records) != expected:
                raise ValueError(f'Incomplete {endpoint}: {len(records)} != {expected}')
            if not records:
                raise ValueError(f'Empty {endpoint}')
            return records
    raise ValueError('Pagination safety limit reached')


def pivot_clinical(records, kind):
    frame = pd.DataFrame(records)
    identity = 'patientId' if kind == 'PATIENT' else 'sampleId'
    required = [identity, 'patientId', 'clinicalAttributeId', 'value']
    if not set(required).issubset(frame.columns) or frame[list(set(required) - {'value'})].isna().any().any():
        raise ValueError(f'Missing identifiers in {kind} clinical records')
    if 'studyId' in frame and not frame.studyId.eq('prad_mskcc').all():
        raise ValueError('Wrong cBioPortal study')
    mapping = frame[[identity, 'patientId']].loc[:, ~pd.Index([identity, 'patientId']).duplicated()].drop_duplicates()
    if mapping[identity].duplicated().any():
        raise ValueError('Sample maps to conflicting patients')
    core = frame[[identity, 'clinicalAttributeId', 'value']]
    duplicates = int(core.duplicated().sum())
    core = core.drop_duplicates()
    if core.duplicated([identity, 'clinicalAttributeId']).any():
        raise ValueError('Conflicting duplicate clinical values')
    wide = core.pivot(index=identity, columns='clinicalAttributeId', values='value').reset_index()
    wide.columns.name = None
    if kind == 'SAMPLE':
        wide = wide.merge(mapping, on='sampleId', validate='one_to_one')
    return wide, duplicates


def decode_dfs(status, raw_time):
    missing_status = pd.isna(status) or str(status).strip() in {'', 'NA', 'N/A', '[Not Available]'}
    missing_time = pd.isna(raw_time) or str(raw_time).strip() in {'', 'NA', 'N/A', '[Not Available]'}
    if not missing_status and status not in DFS_CODES:
        raise ValueError(f'Unknown DFS status: {status!r}')
    event = None if missing_status else DFS_CODES[status]
    months = None if missing_time else float(raw_time)
    if months is not None and (not math.isfinite(months) or months < 0):
        raise ValueError(f'Invalid DFS time: {raw_time!r}')
    reason = 'missing_status' if missing_status else 'missing_time' if missing_time else 'zero_time_requires_review' if months == 0 else 'eligible'
    return event, months, reason


def match_clinical(metadata, patient_records, sample_records, sample_inventory):
    patient, p_dups = pivot_clinical(patient_records, 'PATIENT')
    sample, s_dups = pivot_clinical(sample_records, 'SAMPLE')
    inventory = pd.DataFrame(sample_inventory)
    if not {'sampleId', 'patientId'}.issubset(inventory.columns):
        raise ValueError('Missing cBio sample inventory identifiers')
    inventory = inventory[['sampleId', 'patientId']]
    if inventory.isna().any().any() or inventory.sampleId.duplicated().any():
        raise ValueError('Duplicate/missing cBio sample inventory IDs')
    relation = sample[['sampleId', 'patientId']].merge(inventory, on='sampleId', how='left', suffixes=('', '_inventory'), validate='one_to_one')
    if not relation.patientId.eq(relation.patientId_inventory).all():
        raise ValueError('Clinical patient mapping disagrees with sample inventory')
    audit = metadata.copy()
    audit['selected'] = audit.tumor_type.eq('Primary tumor')
    audit['selection_reason'] = audit.selected.map({True: 'primary_tumor', False: 'not_primary_tumor'})
    audit.loc[audit.tumor_type.isna(), 'selection_reason'] = 'missing_tumor_type'
    primary = audit[audit.selected].copy()
    if primary.empty:
        raise ValueError('No primary tumors')
    sample = sample.rename(columns={c: 'sample__' + c for c in sample.columns if c not in {'sampleId', 'patientId'}})
    patient = patient.rename(columns={c: 'patient__' + c for c in patient.columns if c != 'patientId'})
    matched = primary.merge(sample, left_on='sample_id', right_on='sampleId', how='left', validate='one_to_one', indicator=True)
    if not matched._merge.eq('both').all():
        raise ValueError('Unmatched primary GEO sample IDs: ' + ','.join(matched.loc[matched._merge.ne('both'), 'sample_id']))
    matched = matched.drop(columns='_merge').merge(patient, on='patientId', how='left', validate='many_to_one', indicator=True)
    if not matched._merge.eq('both').all():
        raise ValueError('Primary sample without patient clinical record')
    matched = matched.drop(columns='_merge')
    if matched.patientId.duplicated().any():
        raise ValueError('Multiple primary profiles per patient: explicit outcome-independent selection required')
    if 'sample__SAMPLE_TYPE' in matched and not matched['sample__SAMPLE_TYPE'].eq('Primary').all():
        raise ValueError('GEO/cBio primary tumor labels disagree')
    required = ['patient__DFS_STATUS', 'patient__DFS_MONTHS']
    if not set(required).issubset(matched.columns):
        raise ValueError('Patient DFS fields absent')
    decoded = [decode_dfs(s, t) for s, t in zip(matched[required[0]], matched[required[1]], strict=True)]
    matched['dfs_event'] = pd.array([d[0] for d in decoded], dtype='Int8')
    matched['dfs_months'] = pd.array([d[1] for d in decoded], dtype='Float64')
    matched['dfs_eligibility_reason'] = [d[2] for d in decoded]
    matched['dfs_eligible'] = matched.dfs_eligibility_reason.eq('eligible')
    matched['dfs_status'] = matched[required[0]].map(lambda x: x.split(':', 1)[-1] if pd.notna(x) else None)
    matched['dfs_time_unit'] = 'months'
    matched = matched.rename(columns={'patientId': 'patient_id', 'sampleId': 'cbio_sample_id'})
    audit = audit.merge(inventory.rename(columns={'sampleId': 'cbio_sample_id', 'patientId': 'patient_id'}), left_on='sample_id', right_on='cbio_sample_id', how='left', validate='one_to_one')
    return matched, audit, {'identical_patient_records_collapsed': p_dups, 'identical_sample_records_collapsed': s_dups}
