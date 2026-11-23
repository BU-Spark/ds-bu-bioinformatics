"""Convert the CampLab single-cell atlas module table (Supplementary Table 2) to the module CSV format.

The workbook has one column per module. Row 2 holds the module IDs (L1, L2, ...), row 3 the
category, row 4 the annotation, and rows 5 onward the gene symbols. This writes:

  modules/camplab_modules.csv       gene,module (one row per gene per module; the module-scores input)
  modules/camplab_modules_info.csv  module,category,annotation,n_genes

The workbook is parsed with the standard library only, as cambridge_clinical_source.py does.
Gene symbols are kept as written (whitespace stripped); nothing is renamed or dropped except blanks
and exact duplicates within a module, which are counted and reported.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
REL = '{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
MODULE_ID = re.compile(r'^L\d+$')
ID_ROW, CATEGORY_ROW, ANNOTATION_ROW, FIRST_GENE_ROW = 2, 3, 4, 5


def _column_index(reference: str) -> int:
    letters = re.match(r'[A-Z]+', reference).group(0)
    index = 0
    for letter in letters:
        index = index * 26 + ord(letter) - 64
    return index - 1


def _first_sheet_path(book: zipfile.ZipFile) -> str:
    workbook = ET.fromstring(book.read('xl/workbook.xml'))
    first = workbook.find('s:sheets/s:sheet', NS)
    rels = ET.fromstring(book.read('xl/_rels/workbook.xml.rels'))
    for rel in rels:
        if rel.get('Id') == first.get(REL):
            target = rel.get('Target').lstrip('/')
            return target if target.startswith('xl/') else f'xl/{target}'
    raise ValueError('Workbook has no readable first sheet')


def read_rows(path: Path) -> dict[int, dict[int, str]]:
    """Return {row number: {column index: text}} for the first worksheet."""
    with zipfile.ZipFile(path) as book:
        shared = []
        if 'xl/sharedStrings.xml' in book.namelist():
            for item in ET.fromstring(book.read('xl/sharedStrings.xml')).findall('s:si', NS):
                shared.append(''.join(t.text or '' for t in item.iter(f"{{{NS['s']}}}t")))
        sheet = ET.fromstring(book.read(_first_sheet_path(book)))
    rows: dict[int, dict[int, str]] = {}
    for row in sheet.findall('s:sheetData/s:row', NS):
        cells = {}
        for cell in row.findall('s:c', NS):
            kind, value = cell.get('t'), cell.find('s:v', NS)
            if kind == 'inlineStr':
                text = ''.join(t.text or '' for t in cell.iter(f"{{{NS['s']}}}t"))
            elif value is None:
                continue
            elif kind == 's':
                text = shared[int(value.text)]
            else:
                text = value.text or ''
            text = text.strip()
            if text:
                cells[_column_index(cell.get('r'))] = text
        rows[int(row.get('r'))] = cells
    return rows


def parse_modules(rows: dict[int, dict[int, str]]):
    """Return (info rows, gene rows, report) from the parsed worksheet."""
    ids = rows.get(ID_ROW, {})
    columns = sorted(ids)
    if not columns or any(not MODULE_ID.match(ids[c]) for c in columns):
        raise ValueError(f'Row {ID_ROW} must hold module IDs like L1, L2, ...')
    names = [ids[c] for c in columns]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate module IDs in the header row')
    gene_rows, info, duplicates = [], [], 0
    for column, module in zip(columns, names):
        seen = []
        for number in sorted(n for n in rows if n >= FIRST_GENE_ROW):
            gene = rows[number].get(column)
            if gene is None:
                continue
            if gene in seen:
                duplicates += 1
                continue
            seen.append(gene)
        if not seen:
            raise ValueError(f'Module {module} has no genes')
        gene_rows.extend({'gene': gene, 'module': module} for gene in seen)
        info.append({'module': module, 'category': rows.get(CATEGORY_ROW, {}).get(column, ''),
                     'annotation': rows.get(ANNOTATION_ROW, {}).get(column, ''), 'n_genes': len(seen)})
    sizes = [row['n_genes'] for row in info]
    report = {'modules': len(info), 'gene_rows': len(gene_rows),
              'unique_genes': len({row['gene'] for row in gene_rows}),
              'min_size': min(sizes), 'max_size': max(sizes),
              'modules_under_5_genes': sum(size < 5 for size in sizes),
              'duplicates_within_module_removed': duplicates}
    return info, gene_rows, report


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--input', type=Path, required=True, help='Supp_Table_2_Module_List.xlsx')
    parser.add_argument('--output', type=Path, default=Path('modules/camplab_modules.csv'))
    parser.add_argument('--info', type=Path, default=Path('modules/camplab_modules_info.csv'))
    args = parser.parse_args(argv)
    if not args.input.is_file():
        parser.error(f'input not found: {args.input}')
    info, genes, report = parse_modules(read_rows(args.input))
    write_csv(args.output, genes, ['gene', 'module'])
    write_csv(args.info, info, ['module', 'category', 'annotation', 'n_genes'])
    digest = hashlib.sha256(args.input.read_bytes()).hexdigest()
    print(f'Source {args.input.name} SHA-256 {digest}')
    for key, value in report.items():
        print(f'{key}: {value}')
    print(f'Wrote {args.output} and {args.info}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
