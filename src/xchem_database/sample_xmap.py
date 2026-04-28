def sample_xmap(xmap, transform, sample_array):
    xmap.interpolate_values(sample_array, transform)
    return sample_array