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
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
from pyrosetta.rosetta.protocols.relax import FastRelax

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser
import pandas as pd
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from random import choice
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.scoring import ScoreType

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

    :param pose: pose to perform mutation on
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
        new_pose = relax_structure2(perform_mutation(new_pose, pos + 1, aa))
        print(f'{counter}/{len_mutation} mutations inserted')
        counter += 1

    return new_pose


def random_mutation(pose, list_aa_pos: list):
    """Perform a random point mutation at one of the specified residue positions.

    A residue position is randomly selected from the list of allowed mutation
    sites, and a random amino acid from the 20 standard amino acids is chosen
    as the replacement.

    :param pose: Pose to mutate.
    :param list_aa_pos: List of residue positions that are allowed to be mutated.
    :return: A new pose containing the randomly generated mutation.
    """

    # Create a list of all 20 possible amino acids
    amino_acids = ["A","R","N","D","C","Q","E","G","H","I",
                   "L","K","M","F","P","S","T","W","Y","V"]

    # Randomly choose a position
    pos = choice(list_aa_pos)
    current_aa = pose.residue(pos).name1()

    # Choose a replacement that is not the same as the current residue
    choices = [x for x in amino_acids if x != current_aa]
    aa = choice(choices)

    # Perform mutation return pose
    return perform_mutation(pose, pos, aa)


def random_mutation2(pose, list_aa_pos: list):
    """Perform a random point mutation and identify residues within 8 Å
    of the mutated residue.

    :param pose: Pose to mutate.
    :param list_aa_pos: List of residue positions allowed to be mutated.
    :return: Mutated pose and residues within 8 Å of the mutation.
    """

    # Create a list of all 20 possible amino acids
    amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

    # Randomly choose a position
    pos = choice(list_aa_pos)
    current_aa = pose.residue(pos).name1()
    choices = [x for x in amino_acids if x != current_aa]
    aa = choice(choices)

    # Perform mutation
    new_pose = perform_mutation(pose, pos, aa)

    # Select the mutated residue
    target = ResidueIndexSelector(str(pos))

    # Select residues within 20 Å of the mutated residue
    neighbors = NeighborhoodResidueSelector()
    neighbors.set_focus_selector(target)
    neighbors.set_distance(20.0)
    neighbors.set_include_focus_in_subset(True)

    # Get the selected residues
    subset = neighbors.apply(new_pose)

    residues_to_relax = []

    for i in range(1, new_pose.total_residue() + 1):
        if subset[i]:
            residues_to_relax.append(i)

    return new_pose, residues_to_relax



# endregion

def relax_structure(pose_to_relax) -> Pose:
    """Relax a protein structure using Rosetta FastRelax.

    :param pose_to_relax: Pose to be relaxed
    :return: relaxed pose
    """

    # set scoring function
    scorefxn = get_fa_scorefxn()

    # Create a copy so the original pose is left unchanged
    testPose = Pose()
    testPose.assign(pose_to_relax)

    # constrtain backbone to a certain degree
    scorefxn.set_weight(ScoreType.coordinate_constraint, 1.0)

    # Configure MoveMap (Allow interface rigid-body movement & local flexibility)
    movemap = MoveMap()
    movemap.set_bb(True)
    movemap.set_chi(True)
    if testPose.num_jump() > 0:
        movemap.set_jump(
            1, True
        )  # Allows the peptide to shift slightly relative to protein

    # Initialize the FastRelax protocol
    relax = FastRelax()

    # apply move map to relax
    relax.set_movemap(movemap)

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

def fast_minimize_structure(
    pose: Pose, residues_to_relax: list[int]) -> Pose:


    scorefxn = get_fa_scorefxn()
    scorefxn.set_weight(ScoreType.coordinate_constraint, 1.0)
    test_pose = pose.clone()

    # 1. Select target residues and their 8.0A local neighborhood
    target_selector = ResidueIndexSelector()
    for res in residues_to_relax:
        target_selector.append_index(res)

    nbr_selector = NeighborhoodResidueSelector(
        target_selector, 8.0, include_focus_in_subset=True
    )
    nbr_subset = nbr_selector.apply(test_pose)  # Returns a PyRosetta boolean vector

    # 2. Build MoveMap from the selector subset
    movemap = MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(False)
    movemap.set_jump(False)

    for i in range(1, len(nbr_subset) + 1):
        if nbr_subset[i]:
            movemap.set_chi(i, True)
            movemap.set_bb(i, True)

    # 3. Minimize
    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)
    min_mover.min_type("lbfgs_armijo_nonmonotone")
    min_mover.tolerance(0.01)
    min_mover.max_iter(50)

    min_mover.apply(test_pose)

    return test_pose

def fast_minimize_structure2(pose, residues_to_relax, scorefxn):
    test_pose = pose.clone()

    # 1. Select target mutated residue + 8.0A surrounding shell
    target_selector = ResidueIndexSelector()
    for res in residues_to_relax:
        target_selector.append_index(res)

    shell_selector = NeighborhoodResidueSelector(
        target_selector, 8.0, include_focus_in_subset=True
    )

    # 2. STEP A: Local Sidechain Repack (CRITICAL for fixing initial mutation clashes)
    tf = TaskFactory()
    tf.push_back(
        OperateOnResidueSubset(PreventRepackingRLT(), shell_selector, True)
    )
    tf.push_back(
        OperateOnResidueSubset(RestrictToRepackingRLT(), shell_selector, False)
    )

    packer = PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(test_pose)

    # 3. STEP B: Continuous Local Minimization
    nbr_subset = shell_selector.apply(test_pose)
    movemap = MoveMap()
    movemap.set_bb(False)
    movemap.set_chi(False)
    movemap.set_jump(False)

    for i in range(1, len(nbr_subset) + 1):
        if nbr_subset[i]:
            movemap.set_chi(i, True)
            movemap.set_bb(i, True)

    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)
    min_mover.min_type("lbfgs_armijo_nonmonotone")
    min_mover.tolerance(0.01)
    min_mover.max_iter(30)

    min_mover.apply(test_pose)

    return test_pose

def relax_structure2(pose_to_relax: Pose) -> Pose:
    testPose = pose_to_relax.clone()
    scorefxn = get_fa_scorefxn()
    # 2. Configure MoveMap targeting the binding region
    movemap = MoveMap()
    
    # Allow rigid-body movement between LqhIII and NaV1.5
    if testPose.num_jump() > 0:
        movemap.set_jump(1, True)

    # Enable backbone and sidechain flexibility globally OR restricted to interface
    movemap.set_bb(True)
    movemap.set_chi(True)

    # 3. Setup FastRelax
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.set_movemap(movemap)

    # REMOVED: relax.constrain_relax_to_start_coords(True)
    # This allows LqhIII to find its true local energy minimum!

    # 4. Perform Relaxation
    relax.apply(testPose)

    print(f"Bound State Energy: {scorefxn(testPose):.2f} REU")
    return testPose

# region EXPLORATION


def get_residues(pdb_file: str, motif_range=None):
    protein_chain_A = "A"
    protein_chain_B = "B"

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)

    model = structure[0]

    protein_A = model[protein_chain_A]
    protein_B = model[protein_chain_B]

    # Keep only standard amino acids
    protein_A_residues = [r for r in protein_A if r.id[0] == " "]
    protein_B_residues = [r for r in protein_B if r.id[0] == " "]

    # Keep only the motif if specified
    if motif_range is not None:
        start, end = motif_range
        protein_A_residues = [
            r for r in protein_A_residues
            if start <= r.id[1] <= end
        ]

    return protein_A_residues, protein_B_residues


def contacts(pdb_file: str, cutoff: int, motif_range=None):

    protein_A_residues, protein_B_residues = get_residues(pdb_file, motif_range)
    print(f"Protein_A residues : {len(protein_A_residues)}")
    print(f"Protein_B residues   : {len(protein_B_residues)}")



    # COMPUTE MINIMUM HEAVY-CHAIN ATOM DISTANCE BETWEEN EVERY RESIDUE PAIR


    # initialzie matrix of zeros
    contact_matrix = np.zeros((len(protein_A_residues),
                             len(protein_B_residues)))
    
    # for position and amino acid residue in chain A
    for i, A_res in enumerate(protein_A_residues):

        # create a list of atoms that aren't hydrogens
        A_atoms = [a for a in A_res if a.element != "H"]

        # for position and amino acid residue in chain B
        for j, B_res in enumerate(protein_B_residues):

            # create a list of atoms that aren't hydrogens
             B_atoms = [a for a in B_res if a.element != "H"]
            
            # initalize infinity to search for minimum position
             min_dist = np.inf
            
            # for each atom on residue in chain A
             for A_a in A_atoms:
                 # for each atom on the residues in chain B
                 for B_a in B_atoms:
                     
                     # find distance between atoms
                     d = A_a - B_a

                     # find minimm distance
                     if d < min_dist:
                         min_dist = d
            
            # add minimum distance to contact_matrix at position on chain A and 
            # chain B
             contact_matrix[i, j] = min_dist


    # DISPLAY ONLY CONTACTS WITHIN CUTOFF
    plot_matrix = contact_matrix.copy()
    plot_matrix[plot_matrix > cutoff] = np.nan

    return plot_matrix


def contact_map(pdb_file, cutoff, output_name, motif_range=None):


    # find contact region and residues
    plot_matrix = contacts(pdb_file, cutoff, motif_range)
    protein_A_residues, protein_B_residues = get_residues(pdb_file, motif_range)

    # set labels for each protein
    protein_A_labels = [
        f"{r.resname}{r.id[1]}"
        for r in protein_A_residues
    ]

    protein_B_labels = [
        f"{r.resname}{r.id[1]}"
        for r in protein_B_residues
    ]

    # set size of figure
    plt.figure(figsize=(14, 8))

    sns.heatmap(
        plot_matrix,
        cmap="rocket_r",
        vmin=0,
        vmax=cutoff,
        xticklabels=protein_B_labels,
        yticklabels=protein_A_labels,
        linewidths=0.2,
        linecolor="gray",
        cbar_kws={"label":"Minimum heavy-atom distance (Å)"}
    )

    plt.xlabel("LQHIII")
    plt.ylabel("NaV1.5")
    plt.title(f"Protein–Peptide Contacts (< {cutoff} Å)")

    plt.tight_layout()
    plt.savefig(output_name, dpi=300)
    plt.show()



# endregion


def delta_g(pose):

    scorefxn = get_fa_scorefxn()
    iam = InterfaceAnalyzerMover(
    1,          # interface jump
    False,      # tracer output
    scorefxn)
    iam.apply(pose)
    return iam.get_interface_dG()