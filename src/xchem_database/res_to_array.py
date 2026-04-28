import numpy as np

def res_to_array(res):
    poss = []
    atoms = []
    elements = []
    for atom in res.first_conformer():
        pos = atom.pos
        element = atom.element.atomic_number
        # if atom.has_altloc():
        #     raise Exception
        if element == 1:
            continue
        poss.append([pos.x, pos.y, pos.z])
        atoms.append(atom.name)
        elements.append(element)

    return np.array(poss), np.array(atoms), np.array(elements)

