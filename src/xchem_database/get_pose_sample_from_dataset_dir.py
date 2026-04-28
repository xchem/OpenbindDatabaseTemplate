from xchem_database.res_to_array import res_to_array
from xchem_database.tables import known_hit_pose_sample_dtype

import numpy as np



def get_pose_sample_from_dataset_dir(
        res,
        x, y, z,
        idx_pose
):
    centroid = np.array([x, y, z])
    poss, atom, elements = res_to_array(res, )
    com = np.mean(poss, axis=0).reshape((1, 3))
    event_to_lig_com = com - centroid.reshape((1, 3))
    _poss_centered = poss - com
    _rmsd_target = np.copy(_poss_centered) + np.array([22.5, 22.5, 22.5]).reshape(
        (1, 3)) + event_to_lig_com
    size = min(150, _rmsd_target.shape[0])
    atom_array = np.zeros(150, dtype='<U5')
    elements_array = np.zeros(150, dtype=np.int32)
    pose_array = np.zeros((150, 3))
    pose_array[:size, :] = _rmsd_target[:size, :]
    atom_array[:size] = atom[:size]
    elements_array[:size] = elements[:size]

    known_hit_pos_sample = np.array([(
        idx_pose,
        pose_array,
        atom_array,
        elements_array,
    )],
        dtype=known_hit_pose_sample_dtype
    )
    return known_hit_pos_sample