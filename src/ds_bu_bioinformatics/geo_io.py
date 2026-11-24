"""Strict GEO/TSV source readers with bounded numeric-conversion batches."""
import csv
import re

import numpy as np
import pandas as pd

from .expression import open_text

MISSING = {'', 'NA', 'N/A', 'NULL', 'NAN'}


def numeric_block(rows, columns):
    frame = pd.DataFrame(rows, columns=columns)
    for column in columns:
        tokens = frame[column].astype(str).str.strip()
        frame[column] = pd.to_numeric(tokens.mask(tokens.str.upper().isin(MISSING)), errors='raise')
    values = frame.to_numpy(dtype=float)
    if np.isinf(values).any():
        raise ValueError('Infinite source expression')
    return values


def read_geo(path, *, series, platforms, metadata_only=False, collect_features=False, feature_filter=None):
    headers, characteristics, processing, descriptions = {}, [], [], []
    vector_keys = {'!Sample_geo_accession': 'sample_id', '!Sample_title': 'title',
        '!Sample_platform_id': 'platform_id', '!Sample_description': 'source_description',
        '!Sample_source_name_ch1': 'source_name', '!Sample_organism_ch1': 'organism',
        '!Sample_library_strategy': 'library_strategy'}
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if row == ['!series_matrix_table_begin']:
                break
            if not row:
                continue
            if row[0] == '!Sample_characteristics_ch1':
                characteristics.append(row[1:])
            elif row[0] == '!Sample_data_processing':
                processing.append(row[1:])
            elif row[0] == '!Sample_description':
                descriptions.append(row[1:])
            elif row[0] in vector_keys or row[0] == '!Series_geo_accession':
                if row[0] in headers:
                    raise ValueError('Duplicate GEO identity header: ' + row[0])
                headers[row[0]] = row[1:]
        else:
            raise ValueError('GEO matrix start missing')
        if headers.get('!Series_geo_accession') != [series]:
            raise ValueError('Wrong GEO study accession')
        for key in ('!Sample_geo_accession', '!Sample_title', '!Sample_platform_id'):
            if key not in headers:
                raise ValueError('GEO identity header missing: ' + key)
        ids = headers['!Sample_geo_accession']
        n = len(ids)
        if not n or any(not x for x in ids) or len(set(ids)) != n:
            raise ValueError('Empty/duplicate GEO assay IDs')
        if any(len(v) != n for k, v in headers.items() if k.startswith('!Sample_')) or any(len(v) != n for v in [*characteristics, *processing, *descriptions]):
            raise ValueError('GEO metadata length mismatch')
        if not set(headers['!Sample_platform_id']).issubset(platforms) or not processing:
            raise ValueError('Unexpected GEO platform or missing processing provenance')
        records = []
        for i in range(n):
            record = {name: headers[key][i] for key, name in vector_keys.items() if key in headers}
            record['geo_accession'] = record['sample_id']
            for j, values in enumerate(descriptions):
                record['source_description' + (f'_{j + 1}' if j else '')] = values[i]
            for values in characteristics:
                if not values[i].strip():
                    continue
                key, separator, value = values[i].partition(':')
                key = 'geo__' + re.sub('[^a-z0-9]+', '_', key.strip().lower()).strip('_')
                if not separator or key in record:
                    raise ValueError('Malformed/duplicate GEO characteristic')
                record[key] = value.strip()
            records.append(record)
        if next(reader, []) != ['ID_REF', *ids]:
            raise ValueError('GEO expression/sample order mismatch')
        if metadata_only:
            inventory, seen = [], set()
            for row in reader:
                if row == ['!series_matrix_table_end']:
                    break
                if len(row) != n + 1 or not row[0]:
                    raise ValueError('Malformed metadata-source expression table')
                if collect_features:
                    if row[0] in seen:
                        raise ValueError('Duplicate source expression feature')
                    seen.add(row[0])
                    inventory.append(row[0])
            else:
                raise ValueError('Truncated metadata-source GEO table')
            features = pd.DataFrame(index=pd.Index(inventory, name='feature_id')) if collect_features else None
            return pd.DataFrame(records), features, sorted({v for values in processing for v in values})
        features, blocks, batch, seen = [], [], [], set()
        for row in reader:
            if row == ['!series_matrix_table_end']:
                break
            if len(row) != n + 1 or not row[0] or row[0] in seen:
                raise ValueError('Malformed/duplicate GEO expression row')
            seen.add(row[0])
            if feature_filter is not None and row[0] not in feature_filter:
                continue
            features.append(row[0])
            batch.append(row[1:])
            if len(batch) == 4096:
                blocks.append(numeric_block(batch, ids))
                batch = []
        else:
            raise ValueError('Truncated GEO expression table')
        if batch:
            blocks.append(numeric_block(batch, ids))
    values = np.concatenate(blocks) if blocks else np.empty((0, n))
    expression = pd.DataFrame(values, index=pd.Index(features, name='feature_id'), columns=ids)
    return pd.DataFrame(records), expression, sorted({v for values in processing for v in values})


def read_tsv(path, *, feature_filter=None):
    """Unmodified first-column target IDs and numeric sample columns; no implicit joins."""
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        if len(header) < 2 or any(not x for x in header[1:]) or len(set(header[1:])) != len(header[1:]):
            raise ValueError('TSV requires unique nonempty assay columns')
        ids, blocks, batch, seen = [], [], [], set()
        for row in reader:
            if len(row) != len(header) or not row[0] or row[0] in seen:
                raise ValueError('Malformed/duplicate source TSV target')
            seen.add(row[0])
            if feature_filter is not None and row[0] not in feature_filter:
                continue
            ids.append(row[0])
            batch.append(row[1:])
            if len(batch) == 4096:
                blocks.append(numeric_block(batch, header[1:]))
                batch = []
        if batch:
            blocks.append(numeric_block(batch, header[1:]))
    if not blocks:
        raise ValueError('Empty source expression TSV')
    return pd.DataFrame(np.concatenate(blocks), index=pd.Index(ids, name='feature_id'), columns=header[1:])


def tsv_feature_inventory(path):
    """Read original row IDs without allocating a multi-GB numeric matrix."""
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        if len(header) < 2 or len(set(header[1:])) != len(header[1:]) or any(not x for x in header[1:]):
            raise ValueError('TSV requires unique nonempty assay columns')
        features, seen = [], set()
        for row in reader:
            if len(row) != len(header) or not row[0] or row[0] in seen:
                raise ValueError('Malformed/duplicate source TSV target')
            seen.add(row[0]); features.append(row[0])
    return pd.DataFrame(index=pd.Index(features, name='feature_id'))


def align_columns(expression, metadata, aliases):
    """Exact keys only; never guess by partial string or drop an expression assay."""
    candidates = {}
    for r in metadata.to_dict('records'):
        for key in aliases(r):
            if key:
                candidates.setdefault(key, set()).add(r['sample_id'])
    mapping = {}
    for key in expression.columns:
        hits = candidates.get(key, set())
        if len(hits) != 1:
            raise ValueError(f'Expression assay {key!r} has {len(hits)} exact GEO matches')
        mapping[key] = next(iter(hits))
    if len(set(mapping.values())) != len(mapping):
        raise ValueError('Multiple expression columns map to one GEO assay')
    result = expression.rename(columns=mapping)
    return result, pd.DataFrame({'source_expression_id': list(mapping), 'sample_id': list(mapping.values())})
