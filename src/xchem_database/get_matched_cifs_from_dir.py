from xchem_database.get_lig_block_from_path import get_lig_block_from_path
from xchem_database.match_atoms import match_atoms

def get_matched_cifs(cif_paths, known_hit_residue):
    atom_name_array = [atom.name for atom in known_hit_residue.first_conformer() if atom.element.name != 'H']
    matched_paths = []
    for _cif_path in cif_paths:
        block = get_lig_block_from_path(_cif_path)
        match = match_atoms(atom_name_array, block)

        if match:
            matched_paths.append((_cif_path, block, match))
    return matched_paths