from numcodecs import Blosc

comment_dtype = [
    ('idx', '<i4'),
    ('meta_idx', '<i4'),
    ('comment', 'U150'),
]

event_metadata_dtype = [
    ('idx', '<i4'),
    ('event_idx', '<i4'),
    ('res_id', '<U32'),
    ('ligand_data_idx', 'i8'),
    ('pose_data_idx', 'i8'),
    ('system', '<U32'),
    ('dtag', '<U32'),
    ('event_num', 'i8'),
    ('Confidence', '<U32'),
    ('size', '<f4'),
    ('x', '<f4'),
    ('y', '<f4'),
    ('z', '<f4'),
    ('res', '<f4')
]

z_map_sample_dtype = [('idx', '<i4'), ('sample', '<f4', (90, 90, 90))]

xmap_sample_dtype = [('idx', '<i4'), ('sample', '<f4', (90, 90, 90))]

ligand_data_dtype = [
    ('idx', 'i8'),
    ('num_heavy_atoms', 'i8'),
    ('canonical_smiles', '<U300'),
    ('atom_ids', '<U5', (150,)),
    ('connectivity', '?', (150, 150,))
]

known_hit_pose_sample_dtype = [
    ('idx', '<i8'),
    ('positions', '<f4', (150, 3)),
    ('atoms', '<U5', (150,)),
    ('elements', '<i4', (150,)),
    # ('rmsd', '<f4')
]

annotation_dtype = [
    ('idx', '<i4'),
    ('event_map_table_idx', '<i4'),
    ('event_idx', '<i4'),
    ('annotation', '?'),
    ('partition', 'S32')]



ligand_conf_dtype = [
    ('idx', 'i8'),
    ('ligand_data_idx', 'i8'),
    ('num_heavy_atoms', 'i8'),
    ('fragment_canonical_smiles', '<U300'),
    ('ligand_canonical_smiles', '<U300'),
    ('positions', '<f4', (150, 3)),
    ('elements', '<i4', (150,))]
    

valid_smiles_dtype = [
        ('idx', 'i8'),
        ('valid', '?'),
    ]



def _make_event_metadata_table(group):
    table_z_map_sample_metadata = group.create_dataset(
        'z_map_sample_metadata',
        shape=(0,),
        chunks=(1,),
        dtype=event_metadata_dtype
    )
    return table_z_map_sample_metadata


def _make_z_map_sample_table(group):
    table_z_map_sample = group.create_dataset(
        'z_map_sample',
        shape=(0,),
        chunks=(1,),
        dtype=z_map_sample_dtype,
        compressor=Blosc(cname='zstd', clevel=9, shuffle=Blosc.SHUFFLE)
    )

    return table_z_map_sample

def _make_ligand_data_table(group):
    ligand_data_table = group.create_dataset(
        'ligand_data',
        shape=(0,),
        chunks=(1,),
        dtype=ligand_data_dtype
    )
    return ligand_data_table

def _make_known_hit_pose_table(group):
    table_known_hit_pose_sample = group.create_dataset(
        'known_hit_pose',
        shape=(0,),
        chunks=(1,),
        dtype=known_hit_pose_sample_dtype
    )

    return table_known_hit_pose_sample

def _make_annotation_table(group):
    annotation_table = group.create_dataset(
        'annotation',
        shape=(0,),
        chunks=(1,),
        dtype=annotation_dtype
    )

    return annotation_table


def _make_comment_table(group):


    # def _make_z_map_sample_metadata_table(group):
    annotation_table = group.create_dataset(
        'comments',
        shape=(0,),
        chunks=(1,),
        shards=(1000,),
        dtype=comment_dtype
    )
    # return table_z_map_sample_metadata

    return annotation_table


def _make_xmap_sample_table(group):
    table_xmap_sample = group.create_dataset(
        'xmap_sample',
        shape=(0,),
        chunks=(1,),
        dtype=xmap_sample_dtype,
        compressor=Blosc(cname='zstd', clevel=9, shuffle=Blosc.SHUFFLE)
    )
    return table_xmap_sample


def _make_ligand_conf_table(group):
    mol_conf_group = group.create_dataset(
            'ligand_confs',
            shape=(0,),
            chunks=(1,),
            dtype=ligand_conf_dtype
        )
    
def _make_valid_smiles_table(group):
    
    valid_smiles_group = group.create_dataset(
        'valid_smiles',
        shape=(0,),
        chunks=(1,),
        dtype=valid_smiles_dtype
    )