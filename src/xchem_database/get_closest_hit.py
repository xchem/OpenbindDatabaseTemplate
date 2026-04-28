from xchem_database.res_to_array import res_to_array

import numpy as np

def get_closest_hit(centroid, hits):
    distances = {}
    for resid, res in hits.items():
        res_centroid = np.mean(res_to_array(res)[0], axis=0)
        distance = np.linalg.norm(centroid - res_centroid)
        distances[resid] = distance

    closest_resid = min(distances, key=lambda _x: distances[_x])
    return closest_resid, distances[closest_resid]