"""Read the original publication's clinical workbook using cached XLSX values."""
import io
import posixpath
import re
import zipfile
import xml.etree.ElementTree as ET

import pandas as pd

URL = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4588396/supplementaryFiles'
NS = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
HEADERS = ['SubjectID', 'SampleType_1', 'SampleType_2', 'Total Follow up (months)',
           'Biochemical Relapse', 'BCR Type', 'Time to BCR']


def read_clinical_source(path, cohort='cambridge'):
    if cohort not in {'cambridge', 'stockholm'}:
        raise ValueError('Unknown publication cohort')
    with zipfile.ZipFile(path) as archive:
        names = [n for n in archive.namelist() if n.rsplit('/', 1)[-1] == 'mmc2.xlsx']
        if len(names) != 1:
            raise ValueError('Expected exactly one publication mmc2.xlsx')
        with zipfile.ZipFile(io.BytesIO(archive.read(names[0]))) as book:
            strings = []
            if 'xl/sharedStrings.xml' in book.namelist():
                strings = [''.join(n.itertext()) for n in ET.fromstring(book.read('xl/sharedStrings.xml'))]
            workbook = ET.fromstring(book.read('xl/workbook.xml'))
            sheet = next((s for s in workbook.findall('s:sheets/s:sheet', NS)
                          if s.attrib['name'] == 'Suppl. Table 2 clinical data'), None)
            if sheet is None:
                raise ValueError('Publication clinical worksheet absent')
            rid = sheet.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']
            relations = ET.fromstring(book.read('xl/_rels/workbook.xml.rels'))
            target = next(r.attrib['Target'] for r in relations if r.attrib['Id'] == rid)
            target = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)
            xml = ET.fromstring(book.read(target))
            rows = []
            for row in xml.findall('s:sheetData/s:row', NS):
                values = {}
                for cell in row.findall('s:c', NS):
                    col = re.sub(r'\d+', '', cell.attrib['r'])
                    v = cell.find('s:v', NS)
                    value = v.text if v is not None else None
                    if cell.attrib.get('t') == 's' and value is not None:
                        value = strings[int(value)]
                    elif cell.attrib.get('t') == 'inlineStr':
                        value = ''.join(cell.find('s:is', NS).itertext())
                    values[col] = value
                rows.append(values)
    header = next((r for r in rows if r.get('A') == 'SubjectID'), None)
    if header is None or not set(HEADERS).issubset(header.values()):
        raise ValueError('Unexpected publication clinical headers')
    cols = {name: next(c for c, value in header.items() if value == name) for name in HEADERS}
    records = []
    for row in rows:
        subject = str(row.get(cols['SubjectID']) or '')
        selected = subject.startswith('TB') and subject.endswith('_T') if cohort == 'cambridge' else subject.startswith('STKHLM')
        if selected:
            pattern = r'(TB\d{2}\.\d+)_\d+_T' if cohort == 'cambridge' else r'(STKHLM\d+)'
            match = re.fullmatch(pattern, subject)
            if match is None:
                raise ValueError(f'Unexpected {cohort} publication subject ID')
            record = {name: row.get(c) for name, c in cols.items()}
            record['patient_id'] = match[1]
            records.append(record)
    result = pd.DataFrame(records)
    if result.empty:
        raise ValueError('Missing publication primary subject rows')
    collapsed = []
    for patient_id, group in result.groupby('patient_id', sort=False):
        # Multiple tumor regions may share a subject. Collapse only identical
        # clinical records and preserve every original publication row ID.
        if len(group.drop(columns='SubjectID').drop_duplicates()) != 1:
            raise ValueError('Conflicting publication clinical rows for ' + patient_id)
        record = group.iloc[0].to_dict()
        record['SubjectID'] = ';'.join(group.SubjectID)
        collapsed.append(record)
    result = pd.DataFrame(collapsed)
    if not result.SampleType_1.str.strip().str.upper().eq('TUMOUR').all():
        raise ValueError('Non-tumor publication clinical row')
    return result.rename(columns={'SubjectID': 'paper__subject_id', 'SampleType_1': 'paper__sample_type',
        'SampleType_2': 'paper__specimen_type', 'Total Follow up (months)': 'paper__followup_months',
        'Biochemical Relapse': 'paper__bcr_status', 'BCR Type': 'paper__bcr_type', 'Time to BCR': 'paper__bcr_months'})
