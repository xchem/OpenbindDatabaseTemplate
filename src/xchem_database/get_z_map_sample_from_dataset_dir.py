from xchem_database import constants
from xchem_database.load_xmap_from_path import load_xmap_from_path
from xchem_database.sample_xmap import sample_xmap
from xchem_database.tables import z_map_sample_dtype

import numpy as np
import gemmi

def get_z_map_sample_from_dataset_dir(dataset_dir, x, y, z, idx_z_map):
    # Get the zmap
    zmap_path = dataset_dir / constants.PANDDA_ZMAP_TEMPLATE.format(dtag=dataset_dir.name)
    zmap = load_xmap_from_path(zmap_path)

    # Get the transform
    centroid = np.array([x, y, z])
    transform = gemmi.Transform()
    transform.mat.fromlist((np.eye(3) * 0.5).tolist())
    transform.vec.fromlist((centroid - np.array([22.5, 22.5, 22.5])).tolist())

    # Record the 2fofc map sample
    z_map_sample_array = sample_xmap(
        zmap,
        transform,
        np.zeros((90, 90, 90), dtype=np.float32),
    )

    z_map_sample = np.array(
        [(
            idx_z_map,
            z_map_sample_array
        )],
        dtype=z_map_sample_dtype
    )
    return z_map_sample