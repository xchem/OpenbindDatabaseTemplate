import gemmi

def get_structure_from_path(path):
    return gemmi.read_structure(str(path))
