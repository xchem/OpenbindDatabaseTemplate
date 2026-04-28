from xchem_database.get_matched_cifs_from_dir import get_matched_cifs_from_dir
from xchem_database.get_smiles import get_smiles
from xchem_database.get_atom_ids import get_atom_ids
from xchem_database.get_connectivity import get_connectivity
from xchem_database.tables import ligand_data_dtype

import numpy as np

def get_ligand_data_sample_from_dataset_dir(dataset_dir, res, idx_ligand_data):
    compound_dir = dataset_dir / 'ligand_files'

    # Get the associated ligand data
    try:
        matched_cifs = get_matched_cifs_from_dir(
            res,
            compound_dir,

        )
    except:
        return None

    if len(matched_cifs) == 0:
        # rprint(f'NO MATCHED LIGAND DATA!!!!!!')
        return None

    matched_cif = matched_cifs[0]

    smiles = get_smiles(matched_cif)
    # m = Chem.MolFromSmiles(smiles)
    # m2 = m.AddHs()
    atom_ids_array = np.zeros((150,), dtype='<U5')
    atom_ids = get_atom_ids(matched_cif)
    atom_ids_array[:len(atom_ids)] = atom_ids[:len(atom_ids)]
    num_heavy_atoms = len(atom_ids)
    connectivity_array = np.zeros((150, 150), dtype='?')
    connectivity = get_connectivity(matched_cif)
    connectivity_array[:connectivity.shape[0], :connectivity.shape[1]] = connectivity[
                                                                         :connectivity.shape[0],
                                                                         :connectivity.shape[1]]

    ligand_data_sample = np.array([(
        idx_ligand_data,
        num_heavy_atoms,
        smiles,
        atom_ids_array,
        connectivity_array
    )], dtype=ligand_data_dtype)

    return ligand_data_sample