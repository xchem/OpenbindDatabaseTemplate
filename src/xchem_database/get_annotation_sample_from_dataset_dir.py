from xchem_database.get_system_from_dtag import get_system_from_dtag
from xchem_database.tables import annotation_dtype

import numpy as np

def get_annotation_sample_from_dataset_dir(
        dataset_dir, 
        annotation, 
        test_systems, 
        annotation_idx,
        z_map_sample_metadata_idx,
        ):
    # Get the annotation
    if annotation == "High":
        annotation_bool = True
    else:
        annotation_bool = False

    # Get the partition
    system_name = get_system_from_dtag(dataset_dir.name)
    if system_name in test_systems:
        partition = 'test'
    else:
        partition = 'train'

    # Update
    annotation_sample = np.array(
        [
            (
                annotation_idx,
                z_map_sample_metadata_idx,
                z_map_sample_metadata_idx,
                annotation_bool,
                partition
            )
        ],
        dtype=annotation_dtype
    )
    return annotation_sample