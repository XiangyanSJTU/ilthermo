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


def put_unique_ion(name):
    if re.match(r'^(di|tri|mono|bi)(fluoride|potassium|lithium|sodium|sulfate|chloride|tetrafluoroborate|bromide|nitrate|ammonium|cyanimide|heptanoate|butanoate)', name):
        name = re.split(r'^(di|tri|mono|bi)', name)[-1]
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
        'niobium': '[Nb+5]',
        'germanium': '[Ge+2]',
        'calcium': '[Ca+2]',
        'lanthanum': '[La+3]',
        'mercuric': '[Hg+2]',
        'ytterbium': '[Yb+3]',
        'erbium': '[Er+3]',
        'gallium': '[Ga+3]',
        'thallium': '[Tl+]',
        'tetrabromocobaltate(II)': '[Co-2](Br)(Br)(Br)Br',
        'cholinium': 'C[N+](C)(C)CCO',
        'tris(pentafluoroethyl)trifluorophosphate': 'F[P-](F)(F)(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        'trifluoro(perfluoroethyl)phosphate': 'F[P-](F)(F)(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        'trifluorotris(perfluoroethyl)phosphate(V)': 'F[P-](F)(F)(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        'tris(heptafluoropropyl)trifluorophosphate': 'F[P-](F)(F)(C(F)(F)C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F',
        'tris(nonafluorobutyl)trifluorophosphate': 'F[P-](F)(F)(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)',
        'bis(fluorosulfonyl)imide': '[N-](S(=O)(=O)F)S(=O)(=O)F',
        '1,1,2,2,2-pentafluoro-N-[(pentafluoroethyl)sulfonyl]ethanesulfonamide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        '1,1,1-trifluoro-N-[(trifluoromethyl)sulfonyl]methanesulfonamide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis(pentafluoroethyl)sulfonamide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        'ferrocenec': '[Fe+2]',
        'ferrocenea': '[C-]1C=CC=C1',
        'bis(nonafluorobutylsulfonyl)imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
        'heptachlorodialuminate': 'Cl[Al-](Cl)(Cl)[Cl+][Al-](Cl)(Cl)Cl',
        'ammonioacetate': 'NCC(=O)[O-]',
        'bis(1-butyl-3-methylimidazolium)': 'CCCC[n+]1ccn(C)c1',
        'ibuprofenate': 'CC(C)Cc1ccc(C(C)C(=O)[O-])cc1',
        'tricyanomethane': 'C(#N)[C-](C#N)C#N',
        'tetrachloroindate': '[In-](Cl)(Cl)(Cl)Cl',
        'bis(bis(perfluoroethyl)phosphoryl)imide': 'O=P(C(F)(F)C(F)(F)F)([N-]P(=O)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        'hexafluorophosphate(V)': 'F[P-](F)(F)(F)(F)F',
        '1-butyronitrile-3-methylimidazolium': 'N#CCCC[n+]1ccn(C)c1',
        '1-butyronitrile-2,3-dimethylimidazolium': 'N#CCCC[n+]1ccn(C)c(C)1',
        'N-butyronitrilepyridinium': 'N#CCCC[n+]1ccccc1',
        'Butyronitriletrimethylammonium': 'N#CCCC[N+](C)(C)C',
        'trifluorotris(perfluoroethyl)phosphate': 'FC(F)(F)C(F)(F)[P-](F)(F)(F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F',
        '1-butyl-nicotinic acid butyl ester': 'C(CCC)OC(C=1C=[N+](C=CC1)CCCC)=O',
        'terafluoroborate': 'F[B-](F)(F)F',
        'diethanolammonium': 'OCC[NH2+]CCO',
        'N-decyloxymethyl-3-amido-pyridinium': 'C(CCCCCCCCC)OC[N+]1=CC(=CC=C1)C(=O)N',
        'N,N-dimethylformamide': 'C[N+](C)C=O',
        'N-methylpyrrolidone': 'C[N+]1C(CCC1)=O',
        'L-(+)-N-butylleucine ethyl ester': 'CCCC[NH2+][C@@H](CC(C)C)C(=O)OCC',
        'tris(diethylamino)cyclopropenium': 'C(C)N([CH2+]1C(=C1N(CC)CC)N(CC)CC)CC',
        '1-[(1R,2S,5R)-(-)-menthoxymethyl]-3-methylimidazolium': 'C[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '3-ethyl-1-[(1R,2S,5R)-(-)-menthoxymethyl]imidazolium': 'CC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '1-[(1R,2S,5R)-(-)-menthoxymethyl]-3-propylimidazolium': 'CCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '3-butyl-1-[(1R,2S,5R)-(-)-menthoxymethyl]imidazolium': 'CCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '1-[(1R,2S,5R)-(-)-menthoxymethyl]-3-pentylimidazolium': 'CCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '3-hexyl-1-[(1R,2S,5R)-(-)-menthoxymethyl]imidazolium': 'CCCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '1-[(1R,2S,5R)-(-)-menthoxymethyl]-3-octylimidazolium': 'CCCCCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '1-[(1R,2S,5R)-(-)-menthoxymethyl]-3-nonylimidazolium': 'CCCCCCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '3-decyl-1-[(1R,2S,5R)-(-)-menthoxymethyl]imidazolium': 'CCCCCCCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '1-[(1R,2S,5R)-(-)-menthoxymethyl]-3-undecylimidazolium': 'CCCCCCCCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        '3-dodecyl-1-[(1R,2S,5R)-(-)-menthoxymethyl]imidazolium': 'CCCCCCCCCCCC[n+]1ccn(CO[C@@H]2C[C@H](C)CC[C@H]2C(C)C)c1',
        'n-ethyl-4-(n\',n\'-dimethylammonium)pyridinium': 'c1(N(C)C)cc[n+](CC)cc1',
        '1-octyl-3-methylpydridinium': 'C(CCCCCCC)[N+]1=CC(=CC=C1)C',
        '1,1-(((oxybis(ethane-2,1-diyl))bis(oxy))bis(ethane-2,1-diyl))bis(1-methylpyrrolidin-1-ium)': 'O(CCOCC[N+]1(CCCC1)C)CCOCC[N+]1(CCCC1)C',
        '1,1-(3,6,9,12-tetraoxatetradecane-1,14-diyl)bis(1-methylpyrrolidin-1-ium)': 'C(COCCOCCOCCOCC[N+]1(CCCC1)C)[N+]1(CCCC1)C',
        '3,3-(oxybis(ethane-2,1-diyl))bis(1-butyl-1H-imidazol-3-ium)': 'C(C[N+]1=CN(C=C1)CCCC)OCC[N+]1=CN(C=C1)CCCC',
        '3,3-((ethane-1,2-diylbis(oxy))bis(ethane-2,1-diyl))bis(1-butyl-1H-imidazol-3-ium)': 'C(COCC[N+]1=CN(C=C1)CCCC)OCC[N+]1=CN(C=C1)CCCC',
        '3,3-(((oxybis(ethane-2,1-diyl))bis(oxy))bis(ethane-2,1-diyl))bis(1-butyl-1H-imidazol-3-ium)': 'C(COCCOCC[N+]1=CN(C=C1)CCCC)OCC[N+]1=CN(C=C1)CCCC',
        '3,3-(3,6,9,12-tetraoxatetradecane-1,14-diyl)bis(1-butyl-1H-imidazol-3-ium)': 'C(COCCOCCOCC[N+]1=CN(C=C1)CCCC)OCC[N+]1=CN(C=C1)CCCC',
        'di[bis(trifluoromethylsulfonyl)imide]': '[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F',
        '1,12-di(N-methylpyrrolidinium)dodecane': 'C1CCC[N+]1(C)CCCCCCCCCCCC[N+]1(C)CCCC1',
        'L-phenylalanine ethyl ester': 'C(C)OC([C@@H]([NH3+])CC1=CC=CC=C1)=O',
        'L-leucine methyl ester': 'COC([C@@H]([NH3+])CC(C)C)=O',
        'L-leucine ethyl ester': 'C(C)OC([C@@H]([NH3+])CC(C)C)=O',
        'glycine benzyl ester': 'C(C1=CC=CC=C1)OC(C[NH3+])=O',
        'L-phenylglycine methyl ester': 'COC([C@@H]([NH3+])C1=CC=CC=C1)=O',
        'methyl L-phenylalaninate': '[NH3+][C@@H](CC1=CC=CC=C1)C(=O)OC',
        'L-tryptophan ethyl ester': 'C(C)OC([C@@H]([NH3+])CC1=CNC2=CC=CC=C12)=O',
        'L-phenylalanine tert-butyl ester': 'C(C)(C)(C)OC([C@@H]([NH3+])CC1=CC=CC=C1)=O',
        'L-phenylalanine benzyl ester': 'C(C1=CC=CC=C1)OC([C@@H]([NH3+])CC1=CC=CC=C1)=O',
        'chloroindium': 'Cl[In-](Cl)(Cl)Cl',
        'bis(trifluromethylsulfonyl)imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis{(trifluomethyl)sulfonyl}imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis{(trifluoromethyl)sulfonyl}imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis(trifluoromethanesulfonyl)-imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis(triuoromethylsulfonyl)imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis(trifluoromethylsulfonyl)': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        '(bis(trifluoromethyl)sulfonyl)imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'bis{bis[(trifluoromethyl)sulfonyl]imide}': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F',
        'n-butyl-4-(n\',n\'-dimethylammonium)pyridinium': 'c1(N(C)C)cc[n+](CCCC)cc1',
        '3-methyl-3-propyl-3-azabicyclo[3.2.2]nonanium': 'CCC[N+]1(C)CC2CCC(CC2)C1',
        'propylcholinium': 'CCC[N+](C)(C)CCCO',
        '(1R,2S)-(-)-dimethylephedrinium': 'C[N+](C)(C)[C@@H](C)[C@H](O)c1ccccc1',
        'butyronitriletrimethylammonium': 'C[N+](C)(C)CCCC#N',
        '2-(2-methoxyethoxy) ethylsulfate': 'COCCOCCOS(=O)(=O)[O-]',
        'N-(2-acetyloxy)ethyl-N-methylmorpholinium': 'CC(=O)OCC[N+](C)1CCOCC1',
        'ceftriaxone': 'CO/N=C(\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(CSC3=NC(=O)C([O-])=NN3C)CS[C@H]12)c1csc(N)n1',
        '1,2-bis(dimethyldodecylammonium)ethane': 'CCCCCCCCCCCC[N+](C)(C)CC[N+](C)(C)CCCCCCCCCCCC',
        'tetraphenylphosphorane': 'c1ccc([P+](c2ccccc2)(c2ccccc2)c2ccccc2)cc1',
        'benzenemethanaminium, N-dodecyl-N,N-dimethyl,': 'C(CCCCCCCCCCC)[N+](CC1=CC=CC=C1)(C)C',
        'bromdie': '[Br-]',
        'Benzenemethanaminium, N,N-dimethyl-N-(phenylmethyl)-,': 'C[N+](CC1=CC=CC=C1)(CC1=CC=CC=C1)C',
        'ammonium, trioctylpropyl-,': 'C(CCCCCCC)[N+](CCC)(CCCCCCCC)CCCCCCCC',
        'N,N,N-tri(n-propyl)(4-ethoxy-4-oxobutyl)-1-aminium': 'CCC[N+](CCC)(CCC)CCCC(=O)OCC',
        'N,N,N-triethyl(4-ethoxy-4-oxobutyl)-1-aminium': 'CC[N+](CC)(CC)CCCC(=O)OCC',
        'N,N,N-tri(n-butyl)(4-ethoxy-4-oxobutyl)-1-aminium': 'CCCC[N+](CCCC)(CCCC)CCCC(=O)OCC',
        'N,N,N-tri(n-butyl)(4-ethoxy-4-oxobutyl)-1-phosphonium': 'CCCC[P+](CCCC)(CCCC)CCCC(=O)OCC',
        'cholate': 'C[C@H](CCC(=O)[O-])[C@H]1CC[C@H]2[C@H]3[C@H](C[C@H](O)[C@@]21C)[C@@]1(C)CC[C@@H](O)C[C@H]1C[C@H]3O',
        'warfarin': 'CC(=O)C[C@H](c1ccccc1)c1c([O-])oc2ccccc2c1=O',
        'N-chloro-4-methylbenzenesulfonamide': 'Cl[N-]S(=O)(=O)C1=CC=C(C=C1)C',
        'fluorescein': '[O-]c1ccc2c(Oc3cc([O-])ccc3C24OC(=O)c5ccccc45)c1',
        'trifluoromethansulfate': 'FC(F)(F)OS(=O)(=O)[O-]',
        'N1-(2-ethylhexyl)ethane-1,2-diamine': 'CCCCC(CC)CNCC[NH3+]',
        'trifluoro(perfluoropropyl)borate': 'F[B-](F)(F)(C(F)(F)C(F)(F)C(F)(F)F)',
        'trifluoro(perfluorobutyl)borate': 'F[B-](F)(F)(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)',
        'trifluoro(trifluoromethyl)borate': 'F[B-](F)(F)(C(F)(F)F)',
        'trifluoro(perfluoroethyl)borate': 'F[B-](F)(F)(C(F)(F)C(F)(F)F)',
        '2,4,6-trinitrophenol': '[N+](=O)([O-])C1=C(C(=CC(=C1)[N+](=O)[O-])[N+](=O)[O-])[O-]',
        'cyanocyanamide': 'C(#N)[N-]C#N',
        'docusate': 'CCCCC(CC)COC(=O)CC(C(=O)OCC(CC)CCCC)[S]([O-])(=O)=O',
        'rhodamine': 'CCN(CC)c1ccc2c(c1)[o+]c3cc(ccc3c2c4ccccc4C(O)=O)N(CC)CC',
        'ethanolammonium': 'OCC[NH3+]',
        'hexafluorostibate': 'F[Sb-](F)(F)(F)(F)F',
        'tetraethyammonium': 'C(C)[N+](CC)(CC)CC',
        'L-phenylalanine methyl ester': 'COC([C@@H]([NH3+])CC1=CC=CC=C1)=O',
        '1,1,2,2,2-pentafluoro-N-[(1,1,2,2,2-pentafluoroethyl)sulfonyl]ethanesulfonamide': 'FC(C(F)(F)F)(S(=O)(=O)[N-]S(=O)(=O)C(C(F)(F)F)(F)F)F',
        'c1': 'CN(C)c1cc[n+](COCC(OC[n+]2ccc(N(C)C)cc2)COC[n+]2ccc(N(C)C)cc2)cc1',
        'c2': 'Cn1c[n+](CCCCn2c[n+](C)cc2)cc1',
        'c3': 'Cn1cc[n+](CCCCCCn2cc[n+](C)c2)c1',
        'c4': 'Cn1cc[n+](CCCn2cc[n+](C)c2)c1',
        'bis(dicyanamide)': '[N-](C#N)C#N',
        'cyanimide': '[N-](C#N)C#N',
        'di(dicyanamide)': '[N-](C#N)C#N',
        'tris(dicyanoazanide)': '[N-](C#N)C#N',
        'glycine': 'NCC(=O)[O-]',
        '1-metyl-3-propylimidazolium': 'C[N+]1=CN(C=C1)CCC',
        '(2S)-2-(acetylamino)-5-[(aminoiminomethyl)amino]pentanamide': 'C(C)(=O)N[C@H](C(=O)[NH3+])CCCNC=NN',
        'N,N-dimethylacetamide': 'C[NH+](C(C)=O)C',
        'N,N-dimethylethanamine': 'C[NH+](CC)C',
        'N,N,N\',N\'-tetramethylguanidine': 'CN(C)C(=[NH2+])N(C)C',
        '1-allyl-3-methylimizodalium': 'Cn1cc[n+](CC=C)c1',
        'glycollate': 'OCC(=O)[O-]',
        'dinitramide': '[N-]([N+](=O)[O-])[N+](=O)[O-]',
        'ephedrine': 'C[NH2+][C@@H](C)[C@H](O)c1ccccc1',
        '4-dodecylpyridium': 'C(CCCCCCCCCCC)C1=CC=[NH+]C=C1',
        '4-nonylpyridium': 'C(CCCCCCCC)C1=CC=[NH+]C=C1',
        '3,4-benzo-1-azacyclohexane': 'c1cc2[NH2+]CCCc2cc1',
        '(t-4)-bis(2,4-pentanedionato-O,O\')': 'CC(=O)\\C=C(\\C)[O-]',
        'cefpirome': 'NC=1SC=C(N1)/C(/C(=O)N[C@H]1[C@H]2SCC(=C(N2C1=O)C(=O)O)C[N+]1=C2C(=CC=C1)CCC2)=N/OC',
        'isobutyl L-valinate': '[NH3+][C@@H](C(C)C)C(=O)OCC(C)C',
        'O-isopropyl,L-alanine': 'C(C)(C)OC([C@@H]([NH3+])C)=O',
        'phosphonium, tetraethyl-,': 'C(C)[P+](CC)(CC)CC',
        'pyridinium, 4-(1-hexadecylheptadecyl)-1-methyl-,': 'C(CCCCCCCCCCCCCCC)C(CCCCCCCCCCCCCCCC)C1=CC=[N+](C=C1)C',
        '2,2-dimethylheptane-3,5-dionato': 'CCC(=O)[CH-]C(=O)C(C)(C)C',
        '1-butyl-3-propanenitrile imidazolium': 'CCCCn1cc[n+](CCC#N)c1',
        'dioctysulfosuccinate': 'C(CCCCCCC)OC(CC(S(=O)(=O)[O-])C(=O)OCCCCCCCC)=O',
        'bis(nonafluorobutanesulfonyl)imide': 'O=S(=O)([N-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
        'octadecanoic-acid': 'C(CCCCCCCCCCCCCCCCC)(=O)[O-]',
        'cyclohexanamine': 'C1(CCCCC1)[NH3+]',
        'stearic-acid': 'C(CCCCCCCCCCCCCCCCC)(=O)[O-]',
        'decanoic-acid': 'C(CCCCCCCCC)(=O)[O-]',
        'hexadecanoic-acid': 'C(CCCCCCCCCCCCCCC)(=O)[O-]',
        'palmitic-acid': 'CCCCCCCCCCCCCCCC(=O)[O-]',
        '4,6-dimethyl-N-phenylpyrimidin-2-amine': 'CC1=NC(=NC(=C1)C)[NH2+]C1=CC=CC=C1',
        'heptadecanoic-acid': 'C(CCCCCCCCCCCCCCCC)(=O)[O-]',
        'tetradecanoic-acid': 'CCCCCCCCCCCCCC(=O)[O-]',
        'dodecanoic-acid': 'CCCCCCCCCCCC(=O)[O-]',
        'cyclohexylamine': '[NH3+]C1CCCCC1',
        'tridecanoic-acid': 'C(CCCCCCCCCCCC)(=O)[O-]',
        'diundecanoate': 'C(CCCCCCCCCC)(=O)[O-]',
        'lauric-acid': 'C(CCCCCCCCCCC)(=O)[O-]'
    }
    # special cases
    if name in special_dict:
        smiles = special_dict[name]
    else:
        smiles, inchi = name2smilesinchi(name)
    # use rdkit
    print(smiles)
    rdk_mol = Chem.MolFromSmiles(smiles)
    smiles = mol2smiles(rdk_mol)
    print('put ion:', name, smiles)
    charge = get_charge(smiles=smiles)
    if charge == 0:
        raise Exception('put 0 charged ion')
    ion = UniqueIon(charge=charge, name=name, smiles=smiles)
    ion = add_or_query(ion, 'smiles')
    print('charge:', charge, ';ion_id:', ion.id)
    return ion


def put_ion(name, ion):
    ion_ = Ion(name=name, searched=False, unique_ion=ion)
    add_or_query(ion_, 'name')


def put_unique_molecule(molecule_info, cation, anion):
    if cation is None and anion is None:
        print(molecule_info['name'])
        smiles, inchi = name2smilesinchi(molecule_info['name'])
        if not re.match(r'cannot get SMILES, case \d$', smiles):
            if smiles is not None and '.' in smiles:
                if molecule_info['name'] not in [
                    'potassium sodium (2R,3R)-2,3-dihydroxybutanedioate (1:1:1)',
                    '3-(10,11-dihydro-5H-dibenzo[a,d]cyclohepten-5-ylidine)-N,N-dimethyl-1-propanamine hydrochloride',
                    'acetylferrocene', 'promethazine hydrochloride',
                    '3-(10,11-dihydro-5H-dibenzo[b,f]azepin-5-yl)-N,N-dimethylpropan-1-amine hydrochloride (1:1)',
                    '2-chloro-10-(3-dimethylaminopropyl)phenothiazine hydrochloride',
                    'L-lysine monohydrochloride', 'thiamine hydrochloride',
                    'tris(hydroxymethyl)aminomethane hydrochloride',
                    '5-hydroxy-6-methyl-3,4-pyridinedimethanol hydrochloride', 'guanidine acetate'
                ]:
                    print(smiles)
                    raise Exception('mixture smiles for pure compounds')
            smiles = mol2smiles(Chem.MolFromSmiles(smiles))
    else:
        lcm = np.lcm(cation.charge, anion.charge)
        smiles = '.'.join([anion.smiles] * int(abs(lcm/anion.charge)) + [cation.smiles] * int(lcm/cation.charge))
    molecule = UniqueMolecule(
        name=molecule_info['name'],
        cation=cation,
        anion=anion,
        formula=molecule_info['formula'],
        smiles=smiles
    )
    return add_or_query(molecule, 'smiles')

    # molecule_info['id'] = molecule.id


def put_molecule(molecule_info, molecule):
    molecule_ = Molecule(code=molecule_info['idout'], name=molecule_info['name'], unique_molecule=molecule)
    add_or_query(molecule_, 'name')


def put_paper(paper_info):
    paper = Paper(
        year=paper_info['year'],
        title=paper_info['title'],
        author=paper_info['author'],
        journal=paper_info['journal']
    )
    add_or_query(paper, 'title')
    # paper_info['id'] = paper.id
