import csv
import tempfile
import unittest
import zipfile
from pathlib import Path

from ds_bu_bioinformatics.camplab_modules import main, parse_modules, read_rows

MAIN = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
RELS = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'


def workbook(path: Path, rows: list[list[str]]) -> None:
    """Write a minimal one-sheet xlsx using shared strings, like the CampLab table."""
    strings = sorted({value for row in rows for value in row if value})
    index = {value: i for i, value in enumerate(strings)}
    letters = lambda i: chr(65 + i) if i < 26 else chr(64 + i // 26) + chr(65 + i % 26)
    body = []
    for r, row in enumerate(rows, start=1):
        cells = ''.join(f'<c r="{letters(c)}{r}" t="s"><v>{index[v]}</v></c>' for c, v in enumerate(row) if v)
        body.append(f'<row r="{r}">{cells}</row>')
    with zipfile.ZipFile(path, 'w') as book:
        book.writestr('xl/workbook.xml', f'<workbook xmlns="{MAIN}" xmlns:r="{RELS}"><sheets>'
                      '<sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>')
        book.writestr('xl/_rels/workbook.xml.rels', '<Relationships xmlns="http://schemas.openxmlformats.org/'
                      'package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/>'
                      '</Relationships>')
        book.writestr('xl/sharedStrings.xml', f'<sst xmlns="{MAIN}">' +
                      ''.join(f'<si><t>{s}</t></si>' for s in strings) + '</sst>')
        book.writestr('xl/worksheets/sheet1.xml', f'<worksheet xmlns="{MAIN}"><sheetData>{"".join(body)}'
                      '</sheetData></worksheet>')


TABLE = [['Supplementary Table 2'],
         ['L1', 'L2', 'L3'],
         ['Lineage', 'PCa related Pathway', 'Housekeeping'],
         ['Luminal', 'AR pathway', 'Ribosomal'],
         ['KLK3', 'AR', 'RPL41'],
         ['MSMB', ' KLK2 ', 'RPS27'],
         ['KLK3', 'TMPRSS2', ''],
         ['', 'NKX3-1', '']]


class CampLabModuleTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.xlsx = self.root / 'modules.xlsx'
        workbook(self.xlsx, TABLE)

    def test_wide_table_becomes_gene_module_rows(self):
        info, genes, report = parse_modules(read_rows(self.xlsx))
        self.assertEqual([i['module'] for i in info], ['L1', 'L2', 'L3'])
        self.assertEqual(info[1], {'module': 'L2', 'category': 'PCa related Pathway',
                                   'annotation': 'AR pathway', 'n_genes': 4})
        self.assertEqual([g['gene'] for g in genes if g['module'] == 'L2'], ['AR', 'KLK2', 'TMPRSS2', 'NKX3-1'])
        self.assertEqual(report['duplicates_within_module_removed'], 1)   # second KLK3 in L1
        self.assertEqual(report['modules_under_5_genes'], 3)

    def test_gene_shared_by_modules_is_kept_in_each(self):
        rows = [r[:] for r in TABLE]
        rows[5][2] = 'KLK3'
        workbook(self.xlsx, rows)
        _, genes, _ = parse_modules(read_rows(self.xlsx))
        self.assertEqual(sorted(g['module'] for g in genes if g['gene'] == 'KLK3'), ['L1', 'L3'])

    def test_bad_header_is_rejected(self):
        rows = [r[:] for r in TABLE]
        rows[1][1] = 'Module 2'
        workbook(self.xlsx, rows)
        with self.assertRaises(ValueError):
            parse_modules(read_rows(self.xlsx))

    def test_cli_output_loads_in_module_scores_format(self):
        out, info = self.root / 'm.csv', self.root / 'info.csv'
        self.assertEqual(main(['--input', str(self.xlsx), '--output', str(out), '--info', str(info)]), 0)
        with out.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(list(rows[0]), ['gene', 'module'])
        from ds_bu_bioinformatics.module_scores import load_modules
        try:
            modules = load_modules(out)
        except ImportError:
            self.skipTest('module_scores not importable here')
        self.assertEqual(modules['L2'], ['AR', 'KLK2', 'TMPRSS2', 'NKX3-1'])


if __name__ == '__main__':
    unittest.main()
