from pathlib import Path
import yaml
from rich import print as rprint
import pandas as pd
import numpy as np
import zarr
from rdkit import Chem
from rdkit.Chem import AllChem

from xchem_database import constants
from xchem_database.get_system_from_dtag import get_system_from_dtag
from xchem_database.tables import (
    comment_dtype,
    ligand_conf_dtype,
    valid_smiles_dtype,
    _make_event_metadata_table,
    _make_z_map_sample_table,
    _make_ligand_data_table,
    _make_known_hit_pose_table,
    _make_annotation_table,
    _make_comment_table,
    _make_xmap_sample_table,
    _make_ligand_conf_table,
    _make_valid_smiles_table
)
from xchem_database.get_closest_res_from_dataset_dir import get_closest_res_from_dataset_dir
from xchem_database.get_z_map_sample_from_dataset_dir import get_z_map_sample_from_dataset_dir
from xchem_database.get_pose_sample_from_dataset_dir import get_pose_sample_from_dataset_dir
from xchem_database.get_ligand_data_sample_from_dataset_dir import get_ligand_data_sample_from_dataset_dir
from xchem_database.get_annotation_sample_from_dataset_dir import get_annotation_sample_from_dataset_dir
from xchem_database.get_z_map_metadata_sample_from_dataset_dir import get_z_map_metadata_sample_from_dataset_dir
from xchem_database.get_xmap_sample_from_dataset_dir import get_xmap_sample_from_dataset_dir

def collate_pandda_2_event_data(config_path):
    rprint(f'Running collate_database from config file: {config_path}')
    
    # Load config
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Get the systems to assign to the test partition
    test_systems = config['test']['test_systems']

    # Get the path to create the zarr array in
    zarr_path = config['zarr_path']
    root = zarr.open(zarr_path, mode='a')

    # Delete any previous pandda group in archive
    try:
        del root['pandda_2']
    except:
        rprint(f'No PanDDA 2 group!')
    pandda_2_group = root.create_group('pandda_2')

    # Create 2 new tables in group1
    event_metadata_table = _make_event_metadata_table(pandda_2_group)
    z_map_sample_table = _make_z_map_sample_table(pandda_2_group)
    xmap_sample_table = _make_xmap_sample_table(pandda_2_group)
    ligand_data_table = _make_ligand_data_table(pandda_2_group)
    known_hit_pose_table = _make_known_hit_pose_table(pandda_2_group)
    annotation_table = _make_annotation_table(pandda_2_group)
    comment_table = _make_comment_table(pandda_2_group)

    # Loop over PanDDA directories
    event_metadata_idx = 0
    idx_pose = 0
    idx_ligand_data = 0
    annotation_idx = 0
    for pandda_dir in Path('/dls/data2temp01/labxchem/data/2017/lb18145-17/processing/edanalyzer/output/pandda_new_score/panddas_new_score/').glob('*'):
        # Broken system
        if pandda_dir.name == 'TcCS':
            continue
        # Broken system
        if pandda_dir.name != 'PTP1B':
            continue

        # Skip if not a directory
        if not pandda_dir.is_dir():
            continue
        rprint(f'PanDDA dir is: {pandda_dir}')

        # Skip if no inspect table i.e. no events have been manually reviewed
        pandda_inspect_table = pandda_dir / 'analyses' / 'pandda_inspect_events.csv'
        if not pandda_inspect_table.exists():
            rprint(f'\tNO INSPECT TABLE! SKIPPING!')
            continue

        # Get the inspect table
        inspect_table = pd.read_csv(pandda_inspect_table)

        # Iterate the events contained in the inspect table, creating database records
        for _idx, _row in inspect_table.iterrows():
            
            # Unpack the row information
            dtag, event_idx, bdc, conf, viewed, size = _row['dtag'], _row['event_idx'], _row['1-BDC'], _row[
                constants.PANDDA_INSPECT_HIT_CONDFIDENCE], _row[constants.PANDDA_INSPECT_VIEWED], _row[constants.PANDDA_INSPECT_CLUSTER_SIZE]

            # Assign the system
            system = get_system_from_dtag(dtag)
            rprint(f'\tProcessing event: {dtag} {event_idx} {conf}')

            # Don't create records for unviewed events
            if not viewed:
                rprint('\t\tNot Viewed! Skipping!')
                continue

            # Get the (putative) location of the event
            initial_x, initial_y, initial_z = _row['x'], _row['y'], _row['z']

            # Get the PanDDA 2 dataset dir for the event
            dataset_dir = pandda_dir / 'processed_datasets' / dtag

            # Read the processed dataset yaml and pull out the selected model and events
            processed_dataset_yaml = dataset_dir / 'processed_dataset.yaml'
            with open(processed_dataset_yaml, 'r') as f:
                processed_dataset = yaml.safe_load(f)

            selected_model = processed_dataset['Summary']['Selected Model']

            event_distances = {}
            event_centroids = {}
            for event_num, event in processed_dataset['Models'][selected_model]['Events'].items():
                event_centroid = event['Centroid']
                distance = np.linalg.norm(np.array(event_centroid) - np.array([initial_x, initial_y, initial_z]))
                event_distances[event_num] = distance
                event_centroids[event_num] = event_centroid

            if len(event_distances) > 0:
                closest_event_id = min(event_centroids, key=lambda _event_num: event_distances[_event_num])
                x, y, z = event_centroids[closest_event_id]
                rprint(f'Closest event distance is {event_distances[closest_event_id]}')

            else:
                rprint(
                    f'Could not match high confidence ligand {dtag} {event_idx} to an initial event!\n'
                    f'Check model in {dataset_dir} is appropriate!\n'
                    'SKIPPING!'
                )
                continue

            # Get the corresponding residue
            try:
                resid, res, dist = get_closest_res_from_dataset_dir(
                    dataset_dir,
                    x, y, z
                )
            except Exception as e:
                print(e)
                print('Couldn\'t match res! Skipping!')
                continue
            if (conf == 'High') & (dist > 6.0):
                rprint(
                    f'Could not match high confidence ligand {dtag} {event_idx} to a build!\n'
                    f'Check model in {dataset_dir} is appropriate!\n'
                    'SKIPPING!'
                )
                continue

            # Get the z map sample
            z_map_sample = get_z_map_sample_from_dataset_dir(
                dataset_dir,
                x, y, z,
                event_metadata_idx,
            )
            xmap_sample = get_xmap_sample_from_dataset_dir(
                dataset_dir,
                x, y, z,
                event_metadata_idx,
            )

            # If the pose sample if dataset has been modelled
            if conf == 'High':
                pose_sample = get_pose_sample_from_dataset_dir(
                    res,
                    x, y, z,
                    idx_pose
                )

            # Get the ligand data
            ligand_data_sample = get_ligand_data_sample_from_dataset_dir(
                dataset_dir,
                res,
                idx_ligand_data,
            )
            if not ligand_data_sample:
                rprint(f'\t\tNO LIGAND DATA! SKIPPING!')
                continue

            # Get the annotation data
            annotation_sample = get_annotation_sample_from_dataset_dir(
                dataset_dir,
                conf,
                test_systems,
                annotation_idx,
                event_metadata_idx
            )

            # Get the metadata  
            if conf == 'High':
                tmp_idx_pose = idx_pose
            else:
                tmp_idx_pose = -1
            tmp_idx_ligand_data = idx_ligand_data

            metadata_sample = get_z_map_metadata_sample_from_dataset_dir(
                event_metadata_idx,
                event_idx,
                resid,
                tmp_idx_ligand_data,
                tmp_idx_pose,
                system,
                dtag,
                event_idx,
                conf,
                size,
                x,
                y,
                z,
                _row['high_resolution']
            )

            # Update the tables
            event_metadata_table.append(metadata_sample)
            z_map_sample_table.append(z_map_sample)
            xmap_sample_table.append(xmap_sample)
            ligand_data_table.append(ligand_data_sample)

            comment_sample = np.array(
            [
                (
                    event_metadata_idx,
                    event_metadata_idx,
                    _row['Comment']
                )
            ],
            dtype=comment_dtype
        )
            comment_table.append(comment_sample)

            if conf == 'High':
                known_hit_pose_table.append(pose_sample)
                idx_pose += 1
            idx_ligand_data += 1
            annotation_table.append(annotation_sample)
            event_metadata_idx += 1
            annotation_idx += 1

    # Check their validity (smiles are frequently messed up!)
    valid_smiles_group = _make_valid_smiles_table(pandda_2_group)

    df = pd.DataFrame(
        root['pandda_2']['ligand_data'].get_basic_selection(slice(None), fields=['idx', 'canonical_smiles', ]))

    unique_smiles_series = df['canonical_smiles'].unique()

    smiles_validity = {}
    for idx, smiles in enumerate(unique_smiles_series):

        print(f'{idx}/{len(unique_smiles_series)} : {smiles}')
        try:
            m = Chem.MolFromSmiles(smiles)
            m2 = Chem.AddHs(m)
            cids = AllChem.EmbedMultipleConfs(m2, numConfs=10)
            m3 = Chem.RemoveHs(m2)
            embedding = [_conf.GetPositions() for _conf in m3.GetConformers()][0]

            smiles_validity[smiles] = True
        except Exception as e:
            print(e)
            smiles_validity[smiles] = False

    for _idx, _row in df.iterrows():
        smiles = _row['canonical_smiles']
        if smiles_validity[smiles]:
            valid_smiles_group.append(
                np.array(
                    [(_idx, True)],
                    dtype=valid_smiles_dtype
                )
            )
        else:
            valid_smiles_group.append(
                np.array(
                    [(_idx, False)],
                    dtype=valid_smiles_dtype
                )
            )

    # Generate conformations of ligands

    print(f'Generating ligand confs...')
    ligand_data_table = root['pandda_2']['ligand_data']
    try:
        del root['pandda_2']['ligand_confs']
    except Exception as e:
        print(e)
        _make_ligand_conf_table(pandda_2_group)
    mol_conf_group = root['pandda_2']['ligand_confs']

    mol_conf_idx = 0
    

    high_conf = event_metadata_table.get_mask_selection(event_metadata_table['Confidence'] == ['High'])
    num = len(high_conf)
    for j, _z_map_sample_metadata in enumerate(high_conf):
        print(f'{j} / {num}')
        if _z_map_sample_metadata['Confidence'] != 'High':
            continue
    
        _ligand_data = ligand_data_table[_z_map_sample_metadata['ligand_data_idx']]
    
        m = Chem.MolFromSmiles(_ligand_data['canonical_smiles'])
    
        # Embed the ligand confs (with hydrogens)
        m2 = Chem.AddHs(m)
        cids = AllChem.EmbedMultipleConfs(m2, numConfs=50)
        m3 = Chem.RemoveHs(m2)

        # Go over heavy atom embeddings saving them
        for embedding in [_conf.GetPositions() for _conf in m3.GetConformers()]:
            poss = np.zeros((150, 3))
            poss[:embedding.shape[0], :] = embedding[:, :]
            mol_els = np.array(
                [m2.GetAtomWithIdx(_atom_idx).GetAtomicNum() for _atom_idx in [a.GetIdx() for a in m2.GetAtoms()]])
            els = np.zeros(150)
            els[:len(mol_els)] = mol_els[:]
    
            record = np.array([(
                mol_conf_idx,
                _ligand_data['idx'],
                len(mol_els),
                _ligand_data['canonical_smiles'],
                _ligand_data['canonical_smiles'],
                poss,
                els
            )],
                dtype=ligand_conf_dtype)
    
            #
            mol_conf_group.append(record)
            mol_conf_idx += 1

 