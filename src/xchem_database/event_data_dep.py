import re

import numpy as np
import gemmi

from .database import _get_st_hits, _res_to_array, _get_matched_cifs_from_dir, _get_smiles, \
    _get_atom_ids, _get_connectivity, _get_system_from_dtag
from ..datasets.base import _get_structure_from_path, _load_xmap_from_path, _sample_xmap_and_scale, \
    _load_xmap_from_mtz_path, _sample_xmap
from xchem_database.constants import PANDDA_ZMAP_TEMPLATE, PANDDA_INSPECT_MODEL_DIR, PANDDA_INITIAL_MTZ_TEMPLATE



















def _get_close_distances(known_hit_centroid,
                         experiment_hit_results):
    distances = {}
    for j, res in enumerate(experiment_hit_results):
        if not [x for x in res[0].annotations][0]:
            continue
        centroid = np.array([res[0].x, res[0].y, res[0].z])

        distance = np.linalg.norm(centroid - known_hit_centroid)
        distances[j] = distance
    return distances


def _get_close_events(
        known_hit_centroid,
        experiment_hit_results,
        delta=5.0
):
    distances = _get_close_distances(known_hit_centroid, experiment_hit_results)
    # rprint(distances)
    close_events = []
    for j, dist in distances.items():
        if dist < delta:
            close_events.append(experiment_hit_results[j])

    return close_events


def _get_closest_event(
        known_hit_centroid,
        experiment_hit_results
):
    distances = _get_close_distances(known_hit_centroid, experiment_hit_results)

    return experiment_hit_results[min(distances, key=lambda _j: distances[_j])]






















