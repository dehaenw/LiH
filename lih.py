from sklearn.decomposition import PCA
from rdkit import Chem
from rdkit.Chem import rdDistGeom

nitrile_patt = Chem.MolFromSmarts("[!#1][*]#[*]")
allene_patt = Chem.MolFromSmarts("[*]=[!S;!P]=[*]")
sp2_patt = Chem.MolFromSmarts("[!#1]=[CX3,NX2]-[!#1]")


def get_planarity_of_all_rings(mol):
    """
    detect rings by SSSR, then iterate through rings
    note that this assumes same atom order for comparisons with multiple
    rings. AssignBondOrdersFromTemplate fixes this.
    """
    rings = Chem.GetSSSR(mol) # NOTE THAT THIS ASSUMES SAME ATOM ORDERING
    non_planarities = []
    for ring in rings:
        coords = [mol.GetConformer().GetAtomPosition(idx) for idx in ring]
        pca = PCA(n_components=3)
        t_coords = pca.fit_transform(coords)
        non_planarities += [pca.explained_variance_ratio_[2]*3]
    return non_planarities

def get_collinearity(mol,report_sp2=False):
    """
    detect groups that should be colinear,
    like *C#N; *C#C and C=C=C
    """
    matches = mol.GetSubstructMatches(nitrile_patt)
    matches += mol.GetSubstructMatches(allene_patt)
    if report_sp2:
        matches += mol.GetSubstructMatches(sp2_patt)
    non_collinearities = []
    for match in matches:
        coords = [mol.GetConformer().GetAtomPosition(idx) for idx in match]
        pca = PCA(n_components=3)
        t_coords = pca.fit_transform(coords)
        non_collinearities += [pca.explained_variance_ratio_[1]*2]
    return non_collinearities
    
def check_mols(mols,refs,df,col_threshold=1e-2,p_threshold=1e-3):
    fail_total = []
    for i,m in enumerate(mols):
        fail = False
        if df["Planar"][i]>-1:
            ringps = tuple([1 if p<p_threshold else 0 for p in get_planarity_of_all_rings(m)])
            refrps = tuple([1 if p<p_threshold else 0 for p in get_planarity_of_all_rings(refs[i])])
            if ringps != refrps:
                fail=True
        if df["Collinear"][i]==1:
            cols = get_collinearity(m)
            for col in cols:
                if col>col_threshold:
                    fail=True
        fail_total.append(fail)
    return fail_total
    
def renumber(mol,template):
    try:
        m = Chem.RenumberAtoms(mol, mol.GetSubstructMatch(template))
    except:
        print("remap problem")
        m = mol
    return m
