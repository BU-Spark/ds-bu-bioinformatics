"""Frozen GEO platform annotations without merging expression features."""
import csv

import numpy as np
import pandas as pd


def annotations(path, expression, *, platform, gene_column, gene_pattern,
                symbol_column='Gene Symbol', extra_columns=None, species_column=None):
    """Require an exact platform/schema; quarantine ambiguous and conflicting labels."""
    columns = {'ID': 'feature_id', gene_column: 'gene_id', symbol_column: 'gene_symbol',
               **(extra_columns or {})}
    if species_column:
        columns[species_column] = 'species'
    with path.open(encoding='utf-8') as handle:
        correct = False
        for line in handle:
            correct |= line.rstrip() == f'^PLATFORM = {platform}'
            if line.rstrip() == '!platform_table_begin':
                break
        else:
            raise ValueError(f'{platform}: missing platform table')
        reader = csv.reader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)
        header = next(reader, [])
        if not correct or len(set(header)) != len(header) or not set(columns).issubset(header):
            raise ValueError(f'{platform}: unexpected platform/schema')
        positions = [header.index(c) for c in columns]
        rows = []
        for row in reader:
            if row == ['!platform_table_end']:
                break
            if len(row) != len(header):
                raise ValueError(f'{platform}: malformed annotation row')
            rows.append([row[i].strip() for i in positions])
        else:
            raise ValueError(f'{platform}: truncated platform table')
    frame = pd.DataFrame(rows, columns=list(columns.values()))
    if frame.empty or frame.feature_id.eq('').any() or not frame.feature_id.is_unique:
        raise ValueError(f'{platform}: empty/duplicate feature identifiers')
    if not set(expression.index).issubset(frame.feature_id):
        raise ValueError(f'{platform}: expression feature absent from annotation')
    frame = frame.set_index('feature_id').loc[expression.index].reset_index()
    single = (frame.gene_id.str.fullmatch(gene_pattern)
              & frame.gene_symbol.str.fullmatch(r'[^;/,\s]+')
              & ~frame.gene_symbol.str.upper().isin({'NA', 'N/A', 'NULL', '---', '-'}))
    if species_column:
        single &= frame.species.eq('Homo sapiens')
    pairs = frame.loc[single, ['gene_id', 'gene_symbol']].drop_duplicates()
    bad_ids = set(pairs.loc[pairs.gene_id.duplicated(False), 'gene_id'])
    bad_symbols = set(pairs.loc[pairs.gene_symbol.duplicated(False), 'gene_symbol'])
    conflict = single & (frame.gene_id.isin(bad_ids) | frame.gene_symbol.isin(bad_symbols))
    frame['mapping_status'] = np.select([conflict, single],
        ['gene_symbol_identity_conflict', 'mapped'], default='invalid_or_ambiguous_annotation')
    frame['annotation_eligible'] = frame.mapping_status.eq('mapped')
    return frame
