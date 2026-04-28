from xchem_database.tables import z_map_sample_metadata_dtype

import numpy as np

def get_z_map_metadata_sample_from_dataset_dir(
        idx_z_map,
        event_idx,
        resid,
        ligand_data_idx,
        pose_data_idx,
        system,
        dtag,
        event_num,
        conf,
        size,
x,
y,
                    z,
        res
):
    z_map_sample_metadata = np.array(
        [(
            idx_z_map,
            event_idx,
            resid,
            ligand_data_idx,
            pose_data_idx,
            system,
            dtag,
            event_num,
            conf,
            size,
            x,
            y,
            z,
            res
        )],
        dtype=z_map_sample_metadata_dtype
    )

    return z_map_sample_metadata