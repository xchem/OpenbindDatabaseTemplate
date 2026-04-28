import rich as rprint

def match_atoms(atom_name_array, block):
    atom_id_loop = list(block.find_loop('_chem_comp_atom.atom_id'))
    atom_element_loop = list(block.find_loop('_chem_comp_atom.type_symbol'))

    filtered_atom_id_loop = [_x for _x, _el in zip(atom_id_loop, atom_element_loop) if _el != 'H']

    if len(filtered_atom_id_loop) != len(atom_name_array):
        rprint(f"Different number of atoms! No Match!!")
        return None

    match = {}
    for _j, atom_1_id in enumerate([_x for _x, _el in zip(atom_id_loop, atom_element_loop) if _el != 'H']):
        for _k, atom_2_id in enumerate(atom_name_array):
            if atom_1_id == atom_2_id:
                match[_j] = _k

    if len(match) != len(filtered_atom_id_loop):
        rprint(f"Only partial match {len(match)} / {len(filtered_atom_id_loop)}! Skipping!")
        return None

    else:
        return match