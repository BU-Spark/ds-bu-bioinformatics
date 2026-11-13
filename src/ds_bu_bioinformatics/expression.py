"""Strict parsing of the deposited GPL10264 Series Matrix; no normalization."""
from __future__ import annotations
import csv
import gzip
import re
from pathlib import Path
import numpy as np
import pandas as pd


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if path.suffix == '.gz' else path.open(encoding='utf-8', newline='')


def read_geo(path: Path):
    rows = {}
    characteristics = []
    begin = False
    series = None
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')
        for row in reader:
            if not row:
                continue
            if row[0] == '!series_matrix_table_begin':
                begin = True
                break
            if row[0] == '!Series_geo_accession':
                if series is not None:
                    raise ValueError('Duplicate GEO series accession')
                series = row[1:]
            elif row[0] == '!Sample_characteristics_ch1':
                characteristics.append(row[1:])
            elif row[0] in {'!Sample_geo_accession', '!Sample_title', '!Sample_platform_id', '!Sample_data_processing'}:
                if row[0] in rows:
                    raise ValueError(f'Duplicate GEO header: {row[0]}')
                rows[row[0]] = row[1:]
        if series != ['GSE21034']:
            raise ValueError('Expected GEO series GSE21034')
        required = ['!Sample_geo_accession', '!Sample_title', '!Sample_platform_id', '!Sample_data_processing']
        if not begin or any(key not in rows for key in required):
            raise ValueError('Missing GEO header or table begin')
        accessions = rows['!Sample_geo_accession']
        count = len(accessions)
        if not count or len(set(accessions)) != count or any(not x for x in accessions):
            raise ValueError('Duplicate/empty GEO accessions')
        if any(len(values) != count for values in [*rows.values(), *characteristics]):
            raise ValueError('Inconsistent GEO metadata lengths')
        if set(rows['!Sample_platform_id']) != {'GPL10264'}:
            raise ValueError('Expected expression platform GPL10264')
        metadata = []
        for i in range(count):
            record = {'geo_accession': accessions[i], 'title': rows['!Sample_title'][i], 'platform_id': 'GPL10264'}
            for values in characteristics:
                if not values[i].strip():
                    continue
                if ':' not in values[i]:
                    raise ValueError('Malformed GEO characteristic')
                key, value = values[i].split(':', 1)
                key = re.sub('[^a-z0-9]+', '_', key.lower()).strip('_')
                if key in record:
                    raise ValueError(f'Duplicate characteristic: {key}')
                record[key] = value.strip()
            metadata.append(record)
        meta = pd.DataFrame(metadata)
        if 'sample_id' not in meta or meta.sample_id.eq('').any() or meta.sample_id.duplicated().any():
            raise ValueError('Missing/duplicate GEO sample IDs')
        if 'tumor_type' not in meta:
            raise ValueError('Missing tumor_type')
        header = next(reader, [])
        if header != ['ID_REF', *accessions]:
            raise ValueError('Matrix columns differ from GEO accessions/order')
        features, values = [], []
        end = False
        for row in reader:
            if row == ['!series_matrix_table_end']:
                end = True
                break
            if len(row) != count + 1 or not row[0]:
                raise ValueError('Malformed expression row')
            features.append(row[0])
            values.append(row[1:])
        if not end or not features or len(set(features)) != len(features):
            raise ValueError('Truncated/empty matrix or duplicate feature IDs')
    expression = pd.DataFrame(np.asarray(values, dtype=float), index=pd.Index(features, name='feature_id'), columns=accessions)
    if not np.isfinite(expression.to_numpy()).all():
        raise ValueError('Missing/non-finite expression; no implicit imputation')
    processing = sorted(set(rows['!Sample_data_processing']))
    return meta, expression, processing


def write_csv(frame, path: Path, index=False, float_format='%.17g'):
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = {'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None
    # filename='' avoids embedding staging paths in gzip headers.
    if compression:
        import io
        with path.open('wb') as raw, gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding='utf-8', newline='') as handle:
                frame.to_csv(handle, index=index, lineterminator='\n', float_format=float_format)
    else:
        frame.to_csv(path, index=index, lineterminator='\n', float_format=float_format)
