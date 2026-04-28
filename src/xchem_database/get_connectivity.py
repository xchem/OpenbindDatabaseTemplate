from xchem_database.get_atom_ids import get_atom_ids

import numpy as np

def get_connectivity(matched_cif):
    atom_id_array = get_atom_ids(matched_cif)
    block = matched_cif[1]

    id_to_idx = {}
    for j, atom_id in enumerate(atom_id_array):
        id_to_idx[atom_id] = j
    bond_matrix = np.zeros(
        (
            150,
            150
        ),
        dtype='?')
    bond_1_id_loop = list(block.find_loop('_chem_comp_bond.atom_id_1'))
    bond_2_id_loop = list(block.find_loop('_chem_comp_bond.atom_id_2'))
    for _bond_1_id, _bond_2_id in zip(bond_1_id_loop, bond_2_id_loop):
        if _bond_1_id not in atom_id_array:
            continue
        if _bond_2_id not in atom_id_array:
            continue
        _bond_1_idx, _bond_2_idx = id_to_idx[_bond_1_id], id_to_idx[_bond_2_id]
        bond_matrix[_bond_1_idx, _bond_2_idx] = True
        bond_matrix[_bond_2_idx, _bond_1_idx] = True

    return bond_matrix