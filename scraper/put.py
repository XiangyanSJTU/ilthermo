import time
import sys
import re
import numpy as np
sys.path.append('..')
from scraper.smiles import *
from scraper.models import *


class Log:
    logfile = open('ilscraper-%s.log' % time.strftime('%y%m%d-%H%M%S'), 'w')

    @staticmethod
    def write(*args, **kwargs):
        print(*args, **kwargs, file=sys.stderr)
        print(*args, **kwargs, file=Log.logfile)

    @staticmethod
    def flush():
        Log.logfile.flush()


class SearchFailedError(Exception):
    def __init__(self):
        super().__init__()


class SpecialCaseError(Exception):
    def __init__(self, args=None):
        super().__init__(args)


def add_or_query(row, unique_key):
    filter_dict = {unique_key: row.__dict__[unique_key]}

    result = session.query(row.__class__).filter_by(**filter_dict).first()
    if not result:
        session.add(row)
        session.flush()
        return row
    else:
        return result


def put_prp_table(prp_table):
    prp_table_tuple = [i for i in prp_table.items()]
    prp_table_tuple.sort(key=lambda x: x[0])  # sort by name, alphabet order

    if not session.query(Property).first():  # table exists

        session.bulk_save_objects([
            Property(name=line[0]) for line in prp_table_tuple
        ])
        session.commit()

        return dict([(row[0], i) for row, i in zip(prp_table_tuple, range(len(prp_table_tuple)))])  # name-->idx
    else:
        prps = session.query(Property).all()

        return dict([(prp.name, prp.id) for prp in prps])


def put_ion(name):
    if re.match(r'^(di|tri|mono)(fluoride|potassium|lithium|sodium|sulfate)', name):
        name = re.split(r'^(di|tri|mono)', name)[-1]
    special_dict = {
        'lithium': '[Li+]',
        'sodium': '[Na+]',
        'potassium': '[K+]',
        'Tantalum': '[Ta+5]',
        'cesium': '[Cs+]',
        'rubidium': '[Rb+]',
        'beryllium': '[Be+2]',
        'aluminum': '[Al+3]',
        'magnesium': '[Mg+2]',
        'copper': '[Cu+2]',  # copper monochloride has a bug
        'zinc': '[Zn+2]',
        'ferric': '[Fe+3]',
        'nickel': '[Ni+2]',
        'silver': '[Ag+]',
        'chromium': '[Cr+3]',
        'cholinium': 'C[N+](C)(C)CCO',
        'tris(pentafluoroethyl)trifluorophosphate': 'F[P-](F)(F)(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        'bis(fluorosulfonyl)imide': '[N-](S(=O)(=O)F)S(=O)(=O)F',
        '1,1,2,2,2-pentafluoro-N-[(pentafluoroethyl)sulfonyl]ethanesulfonamide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        '1,1,1-trifluoro-N-[(trifluoromethyl)sulfonyl]methanesulfonamide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'ferrocenec': '[Fe+2]',
        'ferrocenea': '[C-]1C=CC=C1',
        'bis(nonafluorobutylsulfonyl)imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
        'heptachlorodialuminate': 'Cl[Al-](Cl)(Cl)[Cl+][Al-](Cl)(Cl)Cl',
        'ammonioacetate': 'NCC(=O)[O-]',
        'bis(1-butyl-3-methylimidazolium)': 'CCCC[n+]1ccn(C)c1',
        'ibuprofenate': 'CC(C)Cc1ccc(C(C)C(=O)[O-])cc1',
        'tricyanomethane': 'C(#N)[C-](C#N)C#N'
    }
    # special cases
    if name in special_dict:
        smiles = special_dict[name]
    else:
        smiles, inchi = name2smilesinchi(name)
    # use rdkit
    rdk_mol = Chem.MolFromSmiles(smiles)
    smiles = mol2smiles(rdk_mol)
    print('put ion:', name, smiles)
    charge = get_charge(smiles=smiles)
    if charge == 0:
        raise Exception('put 0 charged ion')
    ion = Ion(charge=charge, name=name, smiles=smiles)
    ion = add_or_query(ion, 'smiles')
    print('charge:', charge, ';ion_id:', ion.id)
    return ion


def put_ion_search(name, ion):
    ion_ = IonSearch(name=name, searched=False, ion=ion)
    add_or_query(ion_, 'name')


def put_molecule(molecule_info, cation, anion):
    if cation is None and anion is None:
        smiles, inchi = name2smilesinchi(molecule_info['name'])
        if '.' in smiles:
            if molecule_info['name'] not in [
                'potassium sodium (2R,3R)-2,3-dihydroxybutanedioate (1:1:1)',
                '3-(10,11-dihydro-5H-dibenzo[a,d]cyclohepten-5-ylidine)-N,N-dimethyl-1-propanamine hydrochloride',
                'acetylferrocene'
            ]:
                raise Exception('mixture smiles for pure compounds')
    else:
        lcm = np.lcm(cation.charge, anion.charge)
        smiles = '.'.join([cation.smiles] * int(lcm/cation.charge) + [anion.smiles] * int(abs(lcm/anion.charge)))
    molecule = Molecule(
        code=molecule_info['idout'],
        name=molecule_info['name'],
        cation=cation,
        anion=anion,
        formula=molecule_info['formula'],
        smiles=smiles
    )
    return add_or_query(molecule, 'smiles')

    # molecule_info['id'] = molecule.id


def put_paper(paper_info):
    paper = Paper(
        year=paper_info['year'],
        title=paper_info['title'],
        author=paper_info['author'],
        journal=paper_info['journal']
    )
    add_or_query(paper, 'title')
    # paper_info['id'] = paper.id
