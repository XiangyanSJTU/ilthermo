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
        if name in [
            '1-propan-2-d-ol, 2-methyl-', 'isonicotinic acid hydrazide', '2-butanol (d)',
            'methanesulfonic acid, trifluoro-', 'N,N\'-ethylenebis(salicylideneiminato)diaquachromium(III) chloride',
            'diaquabis(4-methylpyridine)iron(3+)  tris[tetrafluoroborate(1-)]',
            '(N,N-diethylethanamine)(dihydrido)(1-methyl-1H-imidazole-.kappa.N3)boron(1+) bis(trifluoromethylsulfonyl)amide',
            '1-pentanol, 5-phenyl-', '2,2\'-(dodecylimino)bis-ethanol N-oxide',
            '4-hydroxy-2-methyl-N-2-pyridinyl-2H-1,2-benzothiazine-3-carboxamide 1,1-dioxide', 'guanidine acetate',
            '(2-methyloyoxyethyl)dimethylpentyloxyammonium acesulfamate', 'thiophene, tetrahydro-3-methyl-, 1,1-dioxide',
            '2-propanol, 1-(2-butoxy-1-methylethoxy)-', '5-hydroxy-3-methyl-1,2,3-oxadiazolium inner salt'
        ]:
            continue
        special_dict = {
            'monopotassium phosphate': 'potassium dihydrogen phosphate',
            '1H-Imidazolium, 3-butyl-1-methyl-, 1,1,2,3,3,3-hexafluoro-1-propanesulfonate (1:1)': '1-butyl-3-methylimidazolium 1,1,2,3,3,3-hexafluoro-1-propanesulfonate',
            'ferrocene': 'ferrocenec ferrocenea',
            'warfarin sodium': 'sodium warfarin',
            'ammonium chloride (NH4Cl)': 'ammonium chloride',
            'niobium chloride (NbCl5)': 'niobium chloride',
            'nitrogen oxide (NO)': 'nitrogen oxide',
            'difluorogermylene': 'germanium difluoride',
            '1H-Imidazolium, 1,3-bis[(hexyloxy)methyl]-, tetrafluoroborate': '1,3-di(hexyloxy)methyl-imidazolium tetrafluoroborate',
            '(oxybis(ethane-2,1-diyl)) bis(diethylsulfonium)dicyanamide': '(oxybis(ethane-2,1-diyl))bis(diethylsulfonium) dicyanamide',
            'methyl L-phenylalaninate, comp. with bis(trifluoromethylsulfonyl)imide': 'methyl L-phenylalaninate bis(trifluoromethylsulfonyl)imide',
            'N,N-dimethylformamide bis 1,1,1-trifluoro-N-((trifluoromethyl)sulfonyl)methanesulfonamide (1:1)': 'N,N-dimethylformamide bis(trifluoromethylsulfonyl)imide',
            'pyridinium, 1,1\',1\'-[1,2,3-propanetriyltris(oxymethylene)]tris[4-(dimethylamino)-, salt with 1,1,1-trifluoro-N-[(trifluoromethyl)sulfonyl]methanesulfonamide (1:3)': 'c1 1,1,1-trifluoro-N-[(trifluoromethyl)sulfonyl]methanesulfonamide',
            'N-chloro-4-methylbenzenesulfonamide, sodium salt (1:1)': 'sodium N-chloro-4-methylbenzenesulfonamide',
            'fluorescein disodium salt': 'disodium fluorescein',
            '1-butanaminium, N.N.N-tributyl-, salt with 2,4,6-trinitrophenol (1:1)': 'N,N,N-tributyl-1-butanaminium 2,4,6-trinitrophenol',
            '1-propanaminium, N,N,N-tripropyl-, salt with 2,4,6-trinitrophenol (1:1)': 'N,N,N-tripropyl-1-propanaminium 2,4,6-trinitrophenol',
            'benzenemethanaminium, N-dodecyl-N,N-dimethyl, chloride': 'N-dodecyl-N,N-dimethyl,benzenemethanaminium chloride',
            '1-hexanamine, N-hexyl-, hydrochloride': 'N-hexyl-1-hexanaminium chloride',
            'N-octyl-1-octaminium hydrochloride': 'N,N-dioctyl-aminium chloride',
            '1-heptanaminmium, N,N,N-triheptyl-, iodide': 'N,N,N-triheptyl-heptanaminium iodide',
            'rhodamine 6G perchlorate': 'rhodamine perchlorate',
            'beta.-methylcholine chloride': 'beta-methylcholine chloride',
            '2,2-dihydroxydiethylamine': 'di(nitroxyethyl)ammonium nitrate',
            'N,N,N-triethyl-2-methoxyethan-1-aminium trifluoro(perfluoropropyl)borate Chemical Formula: C12H22BF10NO': 'N,N,N-triethyl-2-methoxyethan-1-aminium trifluoro(perfluoropropyl)borate',
            '1H-Imidazolium, 3,3\'-(1,4-butanediyl)bis[1-methyl-, salt with 1,1,2,2,2-pentafluoro-N-[(1,1,2,2,2-pentafluoroethyl)sulfonyl]ethanesulfonamide (1:2)' : 'c2 1,1,2,2,2-pentafluoro-N-[(1,1,2,2,2-pentafluoroethyl)sulfonyl]ethanesulfonamide',
            '1-methyl-3-tetradecylimidazolium acetate [C14MIM][OAC]': '1-methyl-3-tetradecylimidazolium acetate',
            'lead(2+) dipropanoate': 'lead(2+) propanate',
            'ephedrine hydrochloride': 'ephedrine chloride',
            'benzothiazolium, 2-[4-(dimethylamino)phenyl]-3,6-dimethyl-, chloride': '2-[4-(dimethylamino)phenyl]-3,6-dimethyl-benzothiazolium chloride',
            'pyridinium, 4-(1-hexadecylheptadecyl)-1-methyl-, chloride': '4-(1-hexadecylheptadecyl)-1-methyl-pyridinium chloride',
            '(t-4)-bis(2,4-pentanedionato-O,O\')beryllium': 'beryllium (t-4)-bis(2,4-pentanedionato-O,O\')',
            'mercury (II) octoate': 'mercury(II) octanate',
            'cefpirome sulfate': 'cefpirome hydrogen sulfate',
            '1H-imidazolium, 3,3\'-(1,6-hexanediyl)bis[1-methyl-, sulfate (1:2)': 'c3 hydrogen sulfate',
            '1H-Imidazolium, 3,3\'-(1,3-propanediyl)bis[4,5-dihydro-1-methyl-, sulfate (1:2)': 'c4 hydrogen sulfate',
            'L-alanine, 1-methylethyl ester, dodecyl sulfate (1:1)': 'O-isopropyl,L-alanine dodecyl sulfate',
            'bis(2,2-dimethylheptane-3,5-dionato)Cu(II)': 'copper 2,2-dimethylheptane-3,5-dionato',
            'octadecanoic acid cyclohexanamine salt (1:1)': 'cyclohexanamine octadecanoic-acid',
            'hexadecanoic acid cyclohexanamine salt (1:1)': 'cyclohexanamine hexadecanoic-acid',
            'stearic acid, thallium(1+) salt': 'thallium(1+) stearic-acid',
            'decanoic acid, mercuric salt': 'mercuric decanoic-acid',
            'palmitic acid, mercury(2+) salt': 'mercuric palmitic-acid',
            'palmitic acid, thallium(1+) salt': 'thallium(1+) palmitic-acid',
            'decanoic acid cyclohexanamine salt (1:1)': 'cyclohexanamine decanoic-acid',
            'heptadecanoic acid cyclohexanamine salt (1:1)': 'cyclohexanamine heptadecanoic-acid',
            'tetradecanoic acid cyclohexanamine salt(1:1)': 'cyclohexanamine tetradecanoic-acid',
            'dodecanoic acid cyclohexanamine salt (1:1)': 'cyclohexanamine dodecanoic-acid',
            'tridecanoic acid cyclohexylamine salt (1:1)': 'cyclohexylamine tridecanoic-acid',
            'lauric acid, thallium(1+) salt': 'thallium(1+) lauric-acid'
        }
        print(name, name=='octadecanoic acid cyclohexanamine')
        if name in special_dict:
            name = special_dict[name]
        if re.match(r'(|[0-9a-zA-Z()\[\]\',\- ]+) (salt |)\(\d:\d\)', name):
            name = re.split(r' (salt |)\(\d:\d\)', name)[0]
        if name in ['1H-purine-2,6-dione, 3,7-dihydro-3,7-dimethyl-']:
            continue
        s = name.split(' ')
        if s[-1] in ['bis(trifluoromethylsulfonyl)imide', 'bis(trifluoromethanesulfonyl)amide', 'tetrafluoroborate',
                     'bis(trifluoromethylsulfonyl)amide', 'bromide', 'bis(perfluoroethylsulfonyl)imide']:
            m_info['cation'], m_info['anion'] = ' '.join(s[0:-1]), s[-1]
        elif len(s) == 1:
            continue
        elif len(s) == 2:
            if re.match(r'(|[0-9a-zA-Z(),\-]+)(methyl|ethyl|propyl|butyl|vinyl|ethylene|hydrogen|carbonyl|butylene)$', s[0]):
                continue
            if re.match(r'(|di|mon|mono)(oxide|acid|glycol|hydrochloride|oxime|ether|hydrate)$', s[1]):
                continue
            m_info['cation'], m_info['anion'] = s
        elif len(s) == 3:
            if (s[1] in ['hydrogen', 'dihydrogen'] and s[2] in ['phosphate', 'adipate', 'sulfate', 'carbonate']) or \
               (re.match(r'(|[0-9a-zA-Z(),\-]+)(methyl|ethyl|propyl|butyl|hexyl|octyl|dodecyl)$', s[1]) and s[2] in ['sulfate', 'phosphonate', 'phosphate']):
                m_info['cation'], m_info['anion'] = s[0], s[1] + ' ' + s[2]
            # elif re.match(r'(|[0-9a-zA-Z(),\-]+)(methyl|ethyl|propyl|butyl|hexyl|octyl)$', s[0]):
                # continue
            elif re.match(r'(|[0-9a-zA-Z(),\-]+)(methyl|ethyl|propyl|butyl|hexyl|octyl)$', s[0]) and s[1] in ['sulfate', 'ammonium']:
                m_info['cation'], m_info['anion'] = s[0] + s[1], s[2]
            elif s[0] == 'N-butyronitrile' and s[1] == 'pyridinium':
                m_info['cation'], m_info['anion'] = s[0] + s[1], s[2]
            elif re.match(r'(di|tetra|nona)hydrate', s[2]):
                m_info['cation'], m_info['anion'] = s[0], s[1]
            elif s[2] == 'amide':
                m_info['cation'], m_info['anion'] = s[0], s[1] + ' ' + s[2]
            elif re.match(r'(|[0-9a-zA-Z(),\-]+)(methyl|ethyl|propyl|butyl|hexyl|octyl)$', s[0]) and \
                    re.match(r'(|[0-9a-zA-Z(),\-]+)(methyl|ethyl|propyl|butyl|hexyl|octyl|phenyl)$', s[1]) and \
                    re.match(r'(ether|carbonate)', s[2]):
                continue
            elif re.match(r'(|[0-9a-zA-Z(),\-]+)(imidazolium|sodium|choline|imidazol-1-ium)', s[0]):
                m_info['cation'], m_info['anion'] = s[0], s[1] + ' ' + s[2]
            elif re.match(r'(|[0-9a-zA-Z(),\-]+)(imidazolium|sodium|ammonium)', s[1]):
                m_info['cation'], m_info['anion'] = s[0] + ' ' + s[1], s[2]
            elif re.match(r'(iodide|perchlorate)', s[2]):
                m_info['cation'], m_info['anion'] = s[0] + ' ' + s[1], s[2]
            else:
                print(name)
                print(s)
                raise Exception('error input')
        else:
            if name == '(oxybis(ethane-2,1-diyl))bis(diethyl sulfonium) bis((trifluoromethyl)sulfonyl) amide':
                m_info['cation'], m_info['anion'] = '(oxybis(ethane-2,1-diyl))bis(diethyl sulfonium)', 'bis((trifluoromethyl)sulfonyl) amide'
            elif name == '(((oxybis(ethane-2,1-diyl))bis(oxy))bis(ethane-2,1-diyl))bis(diethylsulfonium) bis((trifluoromethyl) sulfonyl) amide':
                m_info['cation'], m_info['anion'] = s[0], s[1] + s[2] + s[3]
            elif name == 'isobutyl L-valinate dodecyl sulfate':
                m_info['cation'], m_info['anion'] = s[0] + ' ' + s[1], s[2] + ' ' + s[3]
    return paper_info, molecule_info, search_result_raw
