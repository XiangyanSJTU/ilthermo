import json
import os
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
        '3-(10,11-dihydro-5H-dibenzo[a,d]cyclohepten-5-ylidine)-N,N-dimethyl-1-propanamine hydrochloride': '[H+].[Cl-].CN(C)CCC=C1c2ccccc2CCc3ccccc13',
        'acetylferrocene': '[Fe+2].[C-]1C(C(=O)C)=CC=C1.[C-]1C=CC=C1',
        '2,2\'-[1,2-phenylenebis(nitrilomethylidyne)]bis-phenol': 'Oc1ccccc1CNc1ccccc1NCc1ccccc1O'
    }
    if name in special_dict:
        smiles = special_dict[name]
        inchi = smiles2inchi(smiles)
        return smiles, inchi
    special_dict = {
        '(.+-.)-.alpha.-aminobutyric acid': 'alpha-aminobutyric acid',
    }
    if name in special_dict:
        name = special_dict[name]
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
        return smiles, inchi
    except:
        return None, None