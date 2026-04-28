from xchem_database.get_st_hits import get_st_hits
from xchem_database.get_closest_hit import get_closest_hit
from xchem_database.get_structure_from_path import get_structure_from_path
from xchem_database.get_most_recent_modelled_structure_from_dataset_dir import get_most_recent_modelled_structure_from_dataset_dir
import numpy as np

def get_closest_res_from_dataset_dir(
        dataset_dir,
        x, y, z
):
    st_path = get_most_recent_modelled_structure_from_dataset_dir(dataset_dir)
    st = get_structure_from_path(st_path)
    hits = get_st_hits(st)
    closest_hit_resid, distance = get_closest_hit(np.array([x, y, z]), hits)

    return closest_hit_resid, hits[closest_hit_resid], distance