import os

from pyrosetta import (init, get_score_function, pose_from_pdb, Pose, get_fa_scorefxn)
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector,
    NeighborhoodResidueSelector,
    NotResidueSelector)
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    InitializeFromCommandline,
    IncludeCurrent,
    NoRepackDisulfides,
    PreventRepackingRLT,
    RestrictToRepackingRLT,
    RestrictAbsentCanonicalAASRLT,
    OperateOnResidueSubset)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax

init(options=[
    '-use_input_sc',
    '-input_ab_scheme', 'AHo_Scheme',
    '-ignore_unrecognized_res',
    '-ignore_zero_occupancy', 'false',
    '-load_PDB_components', 'false',
    '-relax:default_repeats', '2',
    '-no_fconfig',
    '-mute', 'all'  
])

class InsertMutationError(Exception):
    pass

# region MUTATION SETUP
def pack(pose, posi, amino, scorefxn):
    """Mutate a specified residue to a target amino acid and locally repack
    neighboring side chains using the provided Rosetta score function.

    The mutation is restricted to the specified residue, while residues within
    the local neighborhood are allowed to repack. All other residues remain
    fixed. The input pose is modified in place.

    :param pose: Pose to mutate
    :param posi: Position on pose to mutate
    :param amino: amino acid to be mutated in pose
    :param scorefxn: Rosetta scoring function used during packing
    """

    # Select Mutate Position
    mut_posi = ResidueIndexSelector()
    mut_posi.set_index(posi)

    # Select Neighbor Position
    nbr_selector = NeighborhoodResidueSelector()
    nbr_selector.set_distance(8.0) # set neighborhood distance to 8 angstroms
    nbr_selector.set_focus_selector(mut_posi) # Build neighborhood around residue
    nbr_selector.set_include_focus_in_subset(True) # include mutant residue in neighborhood
    
    # Select No Design Area
    not_design = NotResidueSelector(mut_posi)

    # The task factory accepts all the task operations
    tf = TaskFactory()

    # These are TaskOperations which are Instructions for packing

    # Intialize rules stated when you init()
    tf.push_back(InitializeFromCommandline())

    # Includes the current side-chain conformation (rotamer) as one of the candidates
    tf.push_back(IncludeCurrent()) 

    # Prevents residues involved in disulfide bonds (cysteine–cysteine bonds) from being repacked.
    tf.push_back(NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = PreventRepackingRLT()
    #True indicates here that we are flipping the selection.  So that we are turning off everything but the CDR and its neighbors.
    prevent_subset_repacking = OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
    tf.push_back(prevent_subset_repacking)

    # Disable design
    # only change identity of mutant AA
    tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(), not_design))

    # Enable design
    # allows amino acid idnetity to change
    aa_to_design = RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(OperateOnResidueSubset(aa_to_design, mut_posi))
    
    # Create Packer
    packer = PackRotamersMover()
    packer.task_factory(tf)

    # Perform The Move
    if not os.getenv("DEBUG"):
      packer.apply(pose)


def perform_mutation(pose, pos, amino):
    """Wrapper function to perform pack without specifying scoring function

    :param pose: pdb to perform mutation on
    :param pos: position to perform mutation
    :param amino: amino acid to muatate in
    :return: pose of muatated pdb
    """

    # Clone it
    mutant_pose = pose.clone()


    # Create the default Rosetta score function
    scorefxn = get_score_function()

    # Perform packing
    pack(mutant_pose, pos, amino, scorefxn)


    # print useful information
    print()
    print('-' * 50)
    print(f'TASK COMPLETE')
    print(f'Successfully mutated {pose.residue(pos).name()} at position {pos} to {mutant_pose.residue(pos).name()}')
    print(f'Orginal Energy {scorefxn(pose)}; New energy: {scorefxn(mutant_pose)}')
    print('-' * 50)
    print('\n')

    return mutant_pose


def insert_mutation(pdb: str, seq_to_mutate: str, mutation: str) -> Pose:
    """Sequentially insert a mutation. Relaxes structures after each mutation. 

    :param pdb: pdb file containing sequence you wish to mutate
    :param seq_to_mutate: the sequence you want to change
    :param mutation: Mutation you wish to insert
    :return: Pose of pdb file with mutation
    """

    # Find lengths
    len_seq = len(seq_to_mutate)
    len_mutation = len(mutation)

    # raise errror if input is wrong
    if len_seq != len_mutation:
        raise InsertMutationError(
        f"Length mismatch: mutation sequence has {len_mutation} residues, "
        f"but the target sequence has {len_seq} residues.")

    # make pose of pdb
    pose = pose_from_pdb(pdb)

    # create a copy of pose so original is unaltered
    new_pose = Pose()
    new_pose.assign(pose)

    # Grab sequence of pose
    pose_seq = pose.sequence()

    # Find starting index of sequence that needs to be mutated
    start = pose_seq.find(seq_to_mutate)

    # intialize counter
    counter = 1

    # For each position in sequence needed to be mutated
    # Replace with new aa
    for pos, aa in zip(range(start, len(mutation) + start), mutation):
        new_pose = relax_structure(perform_mutation(new_pose, pos + 1, aa))
        print(f'{counter}/{len_mutation} mutations inserted')
        counter += 1

    return new_pose

# endregion

def relax_structure(pose_to_relax):
    """Relax a protein structure using Rosetta FastRelax.

    :param pose_to_relax: Pose to be relaxed
    :param output_name: Path to output pdb of pose
    """

    # Create a copy so the original pose is left unchanged
    testPose = Pose()
    testPose.assign(pose_to_relax)

    # Set up relax parameter
    # Create the full-atom Rosetta score function
    scorefxn = get_fa_scorefxn()

    # Initialize the FastRelax protocol
    relax = FastRelax()

    # Use the full-atom score function during relaxation
    relax.set_scorefxn(scorefxn)

    # Keep the relaxed structure close to the starting coordinates
    relax.constrain_relax_to_start_coords(True)

    # Perform energy minimization and side-chain optimization
    relax.apply(testPose)

    # Print energy after relaxed
    print(f"Relaxed energy: {scorefxn(testPose)}")

    # rename to your desired relaxed structure name

    return testPose
