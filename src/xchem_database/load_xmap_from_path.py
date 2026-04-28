import gemmi

def load_xmap_from_path(path):
    ccp4 = gemmi.read_ccp4_map(str(path))
    ccp4.setup(0.0)
    m = ccp4.grid

    return m