from xchem_database import constants
from xchem_database.load_xmap_from_mtz_path import load_xmap_from_mtz_path
from xchem_database.sample_xmap import sample_xmap
from xchem_database.tables import xmap_sample_dtype

import gemmi
import numpy as np

def get_xmap_sample_from_dataset_dir(dataset_dir, x, y, z, idx_z_map):
    # Get the zmap
    xmap_path = dataset_dir / constants.PANDDA_INITIAL_MTZ_TEMPLATE.format(dtag=dataset_dir.name)
    xmap = load_xmap_from_mtz_path(xmap_path)

    # Get the transform
    centroid = np.array([x, y, z])
    transform = gemmi.Transform()
    transform.mat.fromlist((np.eye(3) * 0.5).tolist())
    transform.vec.fromlist((centroid - np.array([22.5, 22.5, 22.5])).tolist())

    # Record the 2fofc map sample
    xmap_sample_array = sample_xmap(
        xmap,
        transform,
        np.zeros((90, 90, 90), dtype=np.float32),
    )

    xmap_sample = np.array(
        [(
            idx_z_map,
            xmap_sample_array
        )],
        dtype=xmap_sample_dtype
    )
    return xmap_sample
