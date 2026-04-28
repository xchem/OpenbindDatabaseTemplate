from xchem_database.get_fragment_mol_from_dataset_cif_path import get_fragment_mol_from_dataset_cif_path

from rdkit import Chem

def get_smiles(matched_cif):
    mol = get_fragment_mol_from_dataset_cif_path(matched_cif[0])
    smiles = Chem.MolToSmiles(mol)
    return smiles

