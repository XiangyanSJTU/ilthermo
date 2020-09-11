import requests
import sys
import re

sys.path.append('..')
from scraper.put import *


def get_page(url, params={}, try_times=5):
    while try_times > 0:
        try:
            r = requests.get(url, params=params)
            break
        except ConnectionError:
            Log.write('Connection error. trying...')
            try_times -= 1

    if try_times <= 0:
        raise ConnectionAbortedError()
        return None
    else:
        return r.text


def get_prp_table(prp_url, try_times=5):
    """ Get property --> prpcode table.
    """
    try:
        prp_table_raw = get_page(prp_url, try_times=try_times)
    except ConnectionAbortedError:
        Log.write('Cannot get prp table. Using local version instead...')
        prp_table_raw = ''.join([line for line in open('ilprpls.json', 'r').readlines()])  # using local version instead

    prp_table_json = json.loads(prp_table_raw)
    prp_table = {}

    for plist in prp_table_json['plist']:
        prp_table.update(dict(zip(plist['name'], plist['key'])))

    return prp_table


def get_paper_table(search_url, params, try_times=5):
    """ Return formatted paper table.
    """
    try:
        search_result_raw = get_page(search_url, params, try_times)
    except ConnectionAbortedError:
        Log.write('Search failed')
        raise SearchFailedError()

    search_result_json = json.loads(search_result_raw)

    try:
        data_header = search_result_json['header']
    except KeyError:
        Log.write('No result:', search_url)
        raise SearchFailedError()

    code_idx = data_header.index('setid')
    ref_idx = data_header.index('ref')
    prp_idx = data_header.index('prp')
    phase_idx = data_header.index('phases')
    cmp1_idx = data_header.index('cmp1')
    cmp1name_idx = data_header.index('nm1')
    cmp2_idx = data_header.index('cmp2')
    cmp2name_idx = data_header.index('nm2')
    cmp3_idx = data_header.index('cmp3')
    cmp3name_idx = data_header.index('nm3')

    paper_table = []

    for line in search_result_json['res']:
        paper_table.append({
            'code': line[code_idx],
            'ref': line[ref_idx],
            'molecule1_code': line[cmp1_idx],
            'molecule1': line[cmp1name_idx],
            'molecule2_code': line[cmp2_idx],
            'molecule2': line[cmp2name_idx] if len(line) > 9 else None,
            'molecule3_code': line[cmp3_idx],
            'molecule3': line[cmp3name_idx] if len(line) > 10 else None,
            'property': line[prp_idx],
            'phase': line[phase_idx],
        })

    return paper_table


def get_paper_molecule_info(search_url, params, try_times=5):
    parse_success = False
    while not parse_success:
        try:
            search_result_raw = get_page(search_url, params, try_times)
        except ConnectionAbortedError:
            Log.write('Get data failed')
            raise SearchFailedError()

        try:
            search_data = json.loads(search_result_raw)
            parse_success = True
        except json.JSONDecodeError:
            Log.write('Cannot parse. Try downloading again...')

    s = re.split(r' \(|\) ', search_data['ref']['full'])
    paper_info = {
        'title': search_data['ref']['title'],
        'author': s[0],
        'year': int(s[1]),
        'journal': s[2]
    }

    molecule_info = search_data['components']
    for m_info in molecule_info:
        m_info['formula'] = re.sub('\</?SUB\>', '', m_info['formula'])

    for m_info in molecule_info:
        name = m_info['name']
        special_dict = {
            'monopotassium phosphate': 'potassium dihydrogen phosphate',
            '1H-Imidazolium, 3-butyl-1-methyl-, 1,1,2,3,3,3-hexafluoro-1-propanesulfonate (1:1)': '1-butyl-3-methylimidazolium 1,1,2,3,3,3-hexafluoro-1-propanesulfonate',
            'ferrocene': 'ferrocenec ferrocenea',
        }
        if name in special_dict:
            name = special_dict[name]
        s = name.split(' ')
        if len(s) == 1:
            '''
            or m_info['name'] in [
            'carbon dioxide',
            'dimethyl sulfoxide',
            'ethyl acetate',
            'dibutyl ether'
            ]:
            '''
            continue
        elif len(s) == 2:
            if re.match(r'(|[0-9a-zA-Z\-]+)(methyl|ethyl|propyl|butyl|vinyl|ethylene|hydrogen)$', s[0]):
                continue
            if re.match(r'(|di|mon)(oxide|acid|glycol|hydrochloride)', s[1]):
                continue
            m_info['cation'], m_info['anion'] = s
        elif len(s) == 3:
            if (s[1] in ['hydrogen', 'dihydrogen'] and s[2] in ['phosphate', 'adipate', 'sulfate']) or \
               (re.match(r'(|[0-9a-zA-Z()\-]+)(methyl|ethyl|propyl|butyl|hexyl|octyl)$', s[1]) and s[2] in ['sulfate', 'phosphonate', 'phosphate']):
                m_info['cation'], m_info['anion'] = s[0], s[1] + ' ' + s[2]
            elif name == 'isonicotinic acid hydrazide':
                continue
            else:
                print(s)
                raise Exception('error input')
    return paper_info, molecule_info, search_result_raw
