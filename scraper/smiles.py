import json
import os
import re
from rdkit.Chem import AllChem as Chem
import pybel
from urllib.request import urlopen
from pubchempy import get_compounds, get_properties
from chemspipy import ChemSpider


def get_pybel_smiles(smiles):
    py_mol = pybel.readstring("smi", smiles)
    return py_mol.write('can', opt={'n': None}).strip()


def mol2inchi(rdk_mol):
    return Chem.MolToInchi(rdk_mol)


def mol2smiles(rdk_mol):
    return Chem.MolToSmiles(rdk_mol, canonical=True)


def smiles2inchi(smiles):
    rdk_mol = Chem.MolFromSmiles(smiles)
    return mol2inchi(rdk_mol)


def inchi2smiles(inchi):
    rdk_mol = Chem.MolFromInchi(inchi)
    return mol2smiles(rdk_mol)


def get_charge(inchi=None, smiles=None):
        if inchi is not None:
            rdk_mol = Chem.MolFromInchi(inchi)
            py_mol = pybel.readstring("inchi", inchi)
        elif smiles is not None:
            rdk_mol = Chem.MolFromSmiles(smiles)
            py_mol = pybel.readstring("smi", smiles)
        else:
            raise Exception('inchi or smiles is needed as input')
        if Chem.GetFormalCharge(rdk_mol) == py_mol.charge:
            return py_mol.charge
        else:
            raise Exception('The charge from rdkit and pybel is different')


def name2smilesinchi(name):
    special_dict = {
        'hydrogen': '[H][H]',
        'oxygen': 'O=O',
        'nitrogen': 'N#N',
        'n-butan(ol-d)': 'CCCCO[2H]',
        '3-(10,11-dihydro-5H-dibenzo[a,d]cyclohepten-5-ylidine)-N,N-dimethyl-1-propanamine hydrochloride': '[H+].[Cl-].CN(C)CCC=C1c2ccccc2CCc3ccccc13',
        'acetylferrocene': '[Fe+2].[C-]1C(C(=O)C)=CC=C1.[C-]1C=CC=C1',
        'ferrocenylsulfonyl(trifluoromethylsulfonyl)imide': '[Fe+2].[C-]1C(S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F)=CC=C1.[C-]1C=CC=C1',
        '2,2\'-[1,2-phenylenebis(nitrilomethylidyne)]bis-phenol': 'Oc1ccccc1CNc1ccccc1NCc1ccccc1O',
        'rhodium(1+), [(1,2,5,6-.nu.)-1,5-cyclooctadiene][(2R,2\'R,5R,5\'R)-1,1\'-(1,2-phenylene)bis[2,5-dimethylphospholane-kP]]-, tetrafluoroborate(1-) (1:1)': 'cannot get SMILES, case 1',
        'N,N\'-ethylenebis(salicylideneiminato)diaquachromium(III) chloride': 'cannot get SMILES, case 2',
        'diaquabis(4-methylpyridine)iron(3+)  tris[tetrafluoroborate(1-)]': 'cannot get SMILES, case 3',
        '(N,N-diethylethanamine)(dihydrido)(1-methyl-1H-imidazole-.kappa.N3)boron(1+) bis(trifluoromethylsulfonyl)amide': 'cannot get SMILES, case 4',
        'micoflavin': 'cannot get SMILES, case 5',
        '(2-methyloyoxyethyl)dimethylpentyloxyammonium acesulfamate': 'cannot get SMILES, case 6',
        'rel-(1R,2S)-N-methylephedrine': 'CN(C)[C@@H](C)[C@H](O)c1ccccc1',
        '(S)-(2-methoxycarbonyl)pyrrolidinium': 'COC(=O)[C@H]1[NH2+]CCC1',
        'salnaph': 'Oc1ccccc1/C=N/c1cccc2cccc(/N=C/c3ccccc3O)c12',
        '[bis(salicylidene)ethylenediaminato]oxovanadium': '[V-2](=O)235[N+](=CC1=C(C=CC=C1)O2)CC[N+]3=CC4=CC=CC=C4O5',
        '2,2-(4\',4\'\'-dihydroxy)diphenylpropane': 'c1(O)ccc(C(C)(C)c2ccc(O)cc2)cc1',
        '2,2\'-(dodecylimino)bis-ethanol N-oxide': 'C(CCCCCCCCCCC)N(=O)(CCO)CCO',
        '(+-)-carvedilol': 'COc1ccccc1OCCNCC(O)COc2cccc3[nH]c4ccccc4c23',
        '5-hydroxy-3-methyl-1,2,3-oxadiazolium inner salt': 'OC1=CN([NH2+]O1)C'
    }
    if name in special_dict:
        smiles = special_dict[name]
        print(name, smiles)
        if re.match(r'cannot get SMILES, case \d$', smiles):
            return smiles, None
        inchi = smiles2inchi(smiles)
        return smiles, inchi
    special_dict = {
        '(.+-.)-.alpha.-aminobutyric acid': 'alpha-aminobutyric acid',
    }
    if name in special_dict:
        name = special_dict[name]
    if re.match(r'(|[0-9a-zA-Z()+\[\]\',\- ]+)monohydrate', name):
        name = re.split(r' monohydrate', name)[0]
    try:
        url = 'https://opsin.ch.cam.ac.uk/opsin/' + name.replace(' ', '%20')
        info = json.loads(urlopen(url).read().decode('utf8'))
        return info['smiles'], info['stdinchi']
    except:
        a = 1

    try:
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + name + '/smiles'
        smiles = urlopen(url).read().decode('utf8')
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + name + '/inchi'
        inchi = urlopen(url).read().decode('utf8')
        return smiles, inchi
    except:
        a = 1

    try:
        cs_token = os.environ['CHEMSPIDER_TOKEN']
        cs = ChemSpider(cs_token)
        MATCHED = False
        pc_d_list = get_properties('IUPACName,IsomericSMILES,CanonicalSMILES,InChI', name, 'name')
        if len(pc_d_list) == 1:
            d = pc_d_list[0]
            smiles = d.get('CanonicalSMILES') or d.get('IsomericSMILES')
            inchi = d.get('InChI')
            MATCHED = True

        if not MATCHED:
            cs_results = cs.search(name)
            cs_results.wait()
            print('ChemSpider: ', name, cs_results)
            if len(cs_results) == 1:
                result = cs_results[0]
                smiles = result.smiles
                inchi = result.inchi
            else:
                result = cs_results[0]
                smiles = result.smiles
                inchi = result.inchi
        return smiles, inchi
    except:
        return None, None