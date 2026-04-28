from xchem_database import constants

import re

def get_most_recent_modelled_structure_from_dataset_dir(dataset_dir):
    model_dir = dataset_dir / constants.PANDDA_INSPECT_MODEL_DIR
    model_paths = {}
    # print(dataset_dir)
    for path in model_dir.glob('*'):
        fitted_model_regex = 'fitted-v([0-9]*).pdb'
        match = re.match(fitted_model_regex, path.name)

        if match:
            model_paths[int(match[1])] = path
        else:
            model_paths[0] = path

    return model_paths[max(model_paths)]