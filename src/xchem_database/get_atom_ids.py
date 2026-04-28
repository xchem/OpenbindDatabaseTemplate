def get_atom_ids(matched_cif):
    atom_ids = list(matched_cif[1].find_loop('_chem_comp_atom.atom_id'))
    atom_types = list(matched_cif[1].find_loop('_chem_comp_atom.type_symbol'))

    return [_x for _x, _y in zip(atom_ids, atom_types) if _y != 'H']
