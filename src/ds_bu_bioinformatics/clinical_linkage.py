"""Explicit sample-linked clinical supplements; reject conflicting source outcomes."""
import json

import numpy as np
import pandas as pd

from .geo_cohort_pipeline import add_outcomes, decode_bcr, number

RECORD_NAME = 'clinical-linkage/records.csv'
PROVENANCE_NAME = 'clinical-linkage/provenance.json'


def linked_decoder(prefix, fallback):
    def decode(record):
        if record.get(f'linked__{prefix}_present') in (True, 'True', 'true', 1):
            event = record.get(f'linked__{prefix}_event')
            event = None if pd.isna(event) else str(int(event))
            return decode_bcr(event, record.get(f'linked__{prefix}_months'), {'0': 0, '1': 1})
        return fallback(record)
    return decode


def apply(clinical, endpoints, records_path, provenance_path, *, accession, source_definitions=None):
    provenance = json.loads(provenance_path.read_text())
    if provenance.get('study_accession') != accession:
        raise ValueError('Clinical supplement accession differs from this cohort')
    if not isinstance(provenance.get('source'), str) or not provenance['source'].strip():
        raise ValueError('Clinical supplement needs a source citation/description')
    definitions = provenance.get('endpoint_definitions', {})
    records = pd.read_csv(records_path, dtype=str, keep_default_na=False)
    required = {'sample_id', 'endpoint', 'event', 'time', 'time_unit'}
    if records.empty or not required.issubset(records):
        raise ValueError('Clinical supplement needs sample_id, endpoint, event, time, time_unit')
    if records.duplicated(['sample_id', 'endpoint']).any():
        raise ValueError('Duplicate sample/endpoint in clinical supplement')
    if not set(records.sample_id).issubset(clinical.sample_id):
        raise ValueError('Clinical supplement contains unknown or unselected sample IDs; verify crosswalk')
    if not set(records.endpoint).issubset(endpoints):
        raise ValueError('Clinical supplement contains undeclared endpoints; do not relabel outcomes')
    result = clinical.copy()
    source_definitions = source_definitions or {}
    indexed = result.set_index('sample_id')
    for prefix, rows in records.groupby('endpoint', sort=False):
        definition = definitions.get(prefix, {})
        if not all(isinstance(definition.get(key), str) and definition[key].strip()
                   for key in ('definition', 'time_origin', 'censoring', 'sample_linkage')):
            raise ValueError(f'{prefix}: document definition, time_origin, censoring and sample_linkage')
        source_definition = source_definitions.get(prefix, {})
        # Exact documented contracts are deliberately conservative: semantic
        # reinterpretation requires a reviewed adapter change, not free text.
        for field in ('definition', 'time_origin', 'censoring'):
            known = source_definition.get(field)
            has_source_values = clinical[f'{prefix}_event' if field == 'definition' else f'{prefix}_months'].notna().any()
            if known:
                if ' '.join(definition[field].split()) != ' '.join(known.split()):
                    raise ValueError(f'{prefix}: incompatible clinical {field}; source contract must match exactly')
            elif has_source_values:
                raise ValueError(f'{prefix}: unresolved source {field}; cannot mix clinical definitions')
        parsed_events, parsed_months = {}, {}
        for r in rows.itertuples():
            if r.event not in {'0', '1'} or r.time_unit not in {'months', 'years'}:
                raise ValueError('Clinical supplement requires event 0/1 and explicit months/years')
            time = number(r.time)
            if time is None:
                raise ValueError('Clinical supplement row lacks an individual duration')
            months = time * (12 if r.time_unit == 'years' else 1)
            if not np.isfinite(months):
                raise ValueError('Clinical supplement duration overflows months')
            old = indexed.loc[r.sample_id]
            if pd.notna(old[f'{prefix}_event']) and int(old[f'{prefix}_event']) != int(r.event):
                raise ValueError(f'{r.sample_id}: supplement event conflicts with deposited {prefix} status')
            old_months = old[f'{prefix}_months']
            if pd.notna(old_months) and not np.isclose(float(old_months), months, rtol=0, atol=1e-8):
                raise ValueError(f'{r.sample_id}: supplement duration conflicts with deposited {prefix} time')
            parsed_events[r.sample_id] = int(r.event)
            parsed_months[r.sample_id] = months
        result[f'linked__{prefix}_present'] = result.sample_id.isin(parsed_events)
        result[f'linked__{prefix}_event'] = pd.array(result.sample_id.map(parsed_events), dtype='Int8')
        result[f'linked__{prefix}_months'] = pd.array(result.sample_id.map(parsed_months), dtype='Float64')
        decoder = linked_decoder(prefix, endpoints[prefix])
        endpoints = {**endpoints, prefix: decoder}
        result = add_outcomes(result, decoder, prefix)
        for field in ('event', 'duration'):
            column = f'{prefix}_{field}_source'
            if column not in result:
                result[column] = result.get(f'{prefix}_clinical_source', 'deposited source')
            result.loc[result[f'linked__{prefix}_present'], column] = provenance['source']
        result.loc[result[f'linked__{prefix}_present'], f'{prefix}_clinical_source'] = provenance['source']
    return result, endpoints, records, provenance
