"""Utility package to perform engineering of LQHIII

Includes functions to mutate NaV1.5, view contacts, and mutate LQHIII to 
increase binding.


Author: Jimmy Capecci
"""

# region IMPORT LIBRARIES & Initialize PyRosetta

# Import pyrosettta functions and objects
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
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.scoring import ScoreType

# Import other libaries 
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import PDBParser
import pandas as pd
from random import choice, random
import tempfile
from math import exp



# Initialize pyrosetta
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

# endregion


# region CLASSES

class InsertMutationError(Exception):
    pass


class InteractionError(Exception):
    pass


class SaltBridgeError(Exception):
    pass


class HydrogenBondError(Exception):
    pass


class AffinityOptimizerError(Exception):
    pass


class AffinityOptimizer():
    """Optimizes Affinity Between two chains in a pdb file"""


    def __init__(self, pdb: str, scoring_function: str, pos_to_mutate: list, temp: int = 1,
                 cooling_rate: float = 0.95, relax: bool = True, relax_every: int = 1,
                 early_stop_iter: int = 30, number_steps: int = 1000, quiet: bool = False):
        """Initializing affinity optimizer

        :param pdb: Path to the PDB file
        :param scoring_function: Scoring function used to evaluate poses
        :param pos_to_mutate: List of residue positions that can be mutated
        :param temp: Initial temperature indicating how exploratory the run function is, defaults to 1
        :param cooling_rate: How quickly the temperature decreases, defaults to 0.95
        :param relax: Set to True to relax the structure at each step, defaults to True
        :param early_stop_iter: Number of steps without improvement before early stopping the run function, defaults to 30
        :param number_steps: Number of steps in the run function, defaults to 1000
        """

        # initialize all input variable
        self.pdb = pdb
        self.current_pose = pose_from_pdb(pdb)
        self.pos_to_mutate = pos_to_mutate
        self.scoring_function = scoring_function
        self.temp = temp
        self.cooling_rate = cooling_rate
        self.relax = relax
        self.relax_every = max(1, relax_every)
        self.max_no_improve = early_stop_iter
        self.number_steps = number_steps
        self.quiet = quiet
        
        # initialize other important variables that will be tracked
        self.scorefxn = get_fa_scorefxn()
        self.iam = None
        self.best_pose = Pose()
        self.best_pose.assign(self.current_pose)
        self.no_improve_steps = 0
        self.best_score = float('inf')
        self.current_score = None
        self.original_pose = Pose()
        self.original_pose.assign(self.current_pose)

        # if scoring function is selected initialize the following 
        if scoring_function == 'Distance':
            self.func = self.__distance
            self.distances = []
            self.chains = []
            self.positions = []
            self.interactions = []
            self.best_distances = []

        # initialize the following for ΔΔG scoring function
        elif scoring_function == 'DDG':
            self.iam = InterfaceAnalyzerMover(1, False, self.scorefxn)
            self.current_score = self.__interface_dg_score(self.current_pose)
            self.best_score = self.current_score
            self.func = self.__ddg  

        # Raise error if invalid scoring function
        else:
            raise AffinityOptimizerError(f'{scoring_function} is not a valid scoring function. Please'
                                         'enter either "DDG" or "Distance"')


    def __interaction_distance(self, chain1, chain2, pos1, pos2, interaction, new_pose=None):
        """Calculate the shortest potential hydrogen-bond distance between two residues.

        :param pose: PyRosetta Pose containing the residues
        :param chain1: PDB chain ID for the first residue
        :param chain2: PDB chain ID for the second residue
        :param pos1: PDB residue position for the first residue
        :param pos2: PDB residue position for the second residue
        :return: shortest distance between potential hydrogen-bonding atoms in Angstroms
        """
        if new_pose is None:
            new_pose = self.current_pose
        
        # Create a Biopython PDB parser
        parser = PDBParser(QUIET=True)

        structure = self.__get_structure(new_pose)

        # Get the first model from the structure
        model = structure[0]

        # Retrieve the two chains being analyzed
        chainA = model[chain1]
        chainB = model[chain2]

        # Retrieve the residues at the specified PDB residue positions
        res1 = chainA[pos1]
        res2 = chainB[pos2]

        # Get the three-letter amino acid codes for each residue
        aa1 = res1.get_resname()
        aa2 = res2.get_resname()

        # Get the hydrogen-bonding atoms for each residue
        # The returned atom lists correspond to res1 and res2 respectively
        
        if interaction == 'Hydrogen Bond':
            atoms1, atoms2 = self.__hydrogen_bond_atoms(aa1, aa2)
        elif interaction == 'Salt Bridge':
            atoms1, atoms2 = self.__salt_bridge_atoms(aa1, aa2)
        else:
            raise InteractionError('Interaction must be "Hydrogen Bond" or "Salt Bridge"')

        # Store all possible donor/acceptor atom pair distances
        distance = []

        # Compare every possible hydrogen-bonding atom in res1
        # against every possible hydrogen-bonding atom in res2
        for atom1 in atoms1:
            for atom2 in atoms2:

                # Calculate the distance between the two atoms in Angstroms
                distance.append(res1[atom1] - res2[atom2])

        # Return the shortest distance between any possible
        # hydrogen-bonding atom pair
        return float(min(distance))


    def __salt_bridge_atoms(self, aa1, aa2):
        """Identify potential salt-bridge atoms between two amino acids.

        :param aa1: Three-letter amino acid code for the first residue
        :param aa2: Three-letter amino acid code for the second residue
        :return: lists of salt-bridge atoms corresponding to aa1 and aa2
        """

        # Define negatively charged atoms that can participate in salt bridges
        salt_bridge_negative = {
            "ASP": ["OD1", "OD2"],
            "GLU": ["OE1", "OE2"]
        }

        # Define positively charged atoms that can participate in salt bridges
        salt_bridge_positive = {
            "ARG": ["NE", "NH1", "NH2"],
            "LYS": ["NZ"],
            "HIS": ["ND1", "NE2"]
        }

        # Check whether aa1 is positively charged and aa2 is negatively charged
        if (aa1 in salt_bridge_positive and aa2 in salt_bridge_negative):
            return salt_bridge_positive[aa1], salt_bridge_negative[aa2]

        # Check whether aa1 is negatively charged and aa2 is positively charged
        elif (aa1 in salt_bridge_negative and aa2 in salt_bridge_positive):
            return salt_bridge_negative[aa1], salt_bridge_positive[aa2]

        # Raise an error if the residues cannot form a salt bridge
        else:
            raise SaltBridgeError(
                f'{aa1} and {aa2} cannot make a salt bridge'
            )


    def __hydrogen_bond_atoms(self, aa1, aa2):
        """Identify potential hydrogen-bonding atoms between two amino acids.

        :param aa1: Three-letter amino acid code for the first residue
        :param aa2: Three-letter amino acid code for the second residue
        :return: lists of hydrogen-bond donor/acceptor atoms corresponding to aa1 and aa2
        """

        # Define atoms that can donate a hydrogen bond for each amino acid
        hbond_donors = {
            "ARG": ["NE", "NH1", "NH2"],
            "ASN": ["ND2"],
            "CYS": ["SG"],
            "GLN": ["NE2"],
            "HIS": ["ND1", "NE2"],
            "LYS": ["NZ"],
            "SER": ["OG"],
            "THR": ["OG1"],
            "TRP": ["NE1"],
            "TYR": ["OH"]
        }

        # Define atoms that can accept a hydrogen bond for each amino acid
        hbond_acceptors = {
            "ASN": ["OD1"],
            "ASP": ["OD1", "OD2"],
            "CYS": ["SG"],
            "GLN": ["OE1"],
            "GLU": ["OE1", "OE2"],
            "HIS": ["ND1", "NE2"],
            "MET": ["SD"],
            "SER": ["OG"],
            "THR": ["OG1"],
            "TYR": ["OH"]
        }

        # Check whether aa1 can donate and aa2 can accept a hydrogen bond
        if (aa1 in hbond_donors and aa2 in hbond_acceptors):
            return hbond_donors[aa1], hbond_acceptors[aa2]

        # Check whether aa2 can donate and aa1 can accept a hydrogen bond
        elif aa2 in hbond_donors and aa1 in hbond_acceptors:
            return hbond_acceptors[aa1], hbond_donors[aa2]

        # Raise an error if neither residue can form a potential hydrogen bond
        else:
            raise HydrogenBondError(
                f'{aa1} and {aa2} cannot make a hydrogen bond'
            )


    def insert_interaction(self, chain1: str, chain2: str, pos1: int, pos2: int, interaction):
        """inserts interaction into the AffinityOptimizer in order to optimize
           distance.

        :param chain1: Chain for first residues Ex: 'A"
        :param chain2: Chain for second residues Ex: 'B'
        :param pos1: Position in chain for first amino acid
        :param pos2: Position in chain for second amino acid
        :param interaction: specify interaction type Ex: 'Salt Bridge' or 'Hydrogen Bond'
        """

        # Find distance for interaction and append to distances
        distance = self.__interaction_distance(chain1, chain2, pos1, pos2, interaction)
        self.distances.append(distance)

        # keep track of chains, positions, and interactions
        self.chains.append([chain1, chain2])
        self.positions.append([pos1, pos2])
        self.interactions.append(interaction)

        # Initialize best distances and best score for run function
        self.best_distances.append(distance)
        self.best_score = self.__distance_score(self.best_distances)
        self.current_score = self.best_score


    def __distance_score(self, pose_distances: list):
        """finds score of give distances for poses

        :param pose_distances: distances between interactions
        :return: returns score
        """

        # initialize score
        score = 0

        # find score
        for dis in pose_distances:
            error = abs(dis - 2.8)

            score += error ** 2

        return score


    def find_distances(self, pose: Pose =None):
        """Finds distances of all interactions

        :param pose: pose to find distances, defaults to None
        :return: list of distances
        """

        # if no pose inserted use current pose
        if pose is None:
            pose = self.current_pose

        # initalize distances
        dis = []

        # for each interaction find distance for pose
        for i, (dbl_chain, dbl_pos) in enumerate(zip(self.chains, self.positions)):
            dis.append(self.__interaction_distance(dbl_chain[0], dbl_chain[1], dbl_pos[0], dbl_pos[1], self.interactions[i], pose))

        return dis


    def __distance(self, new_pose: Pose):
        """Calculate the distance score for a new pose and evaluate the pose.

        Calculates the distances between the specified interactions in the new
        pose, then uses the distance score and simulated annealing algorithm to
        determine whether the new pose should be accepted.

        :param new_pose: New pose generated from the mutation step
        """

        new_distances = self.find_distances(new_pose)
        new_score = self.__distance_score(new_distances)
        old_score = self.__distance_score(self.distances)
        self.__algorithm(new_pose, new_score, old_score, new_distances)


    def __algorithm(self, new_pose, new_score, old_score, new_distances=None):
        """Apply simulated annealing acceptance and track the best pose.

        Compares the new score to the old score and determines whether to
        accept the new pose. Better-scoring poses are always accepted, while
        worse-scoring poses may be accepted based on the current temperature.
        Updates the current pose, best pose, scores, distances, and temperature
        as needed.

        :param new_pose: New pose generated from the mutation step
        :param new_score: Score of the new pose
        :param old_score: Score of the current pose before mutation
        :param new_distances: Distances calculated for the new pose, defaults to None
        """

        delta_score = new_score - old_score
        accepted = False

        # Standard Metropolis acceptance criterion
        if delta_score < 0:
            accepted = True
        else:
            prob = exp(-delta_score / self.temp)
            if prob > random():
                accepted = True

        # If accepted, update the state to the new pose and score
        if accepted:
            self.current_pose = new_pose
            self.current_score = new_score
            if new_distances is not None:
                self.distances = new_distances

        # Track best score & pose
        if self.current_score < self.best_score:
            self.best_pose = Pose()
            self.best_pose.assign(self.current_pose)
            if self.scoring_function == 'Distance':
                self.best_distances = self.distances.copy()
            self.best_score = self.current_score
            self.no_improve_steps = 0
        else:
            self.no_improve_steps += 1

        # Cool down temperature
        self.temp *= self.cooling_rate


    def __interface_dg_score(self, pose: Pose):
        """Calculate the binding affinity of a protein complex.

        Uses the InterfaceAnalyzerMover to calculate the interface binding
        energy of the given pose.

        :param pose: PyRosetta pose containing the protein complex
        :return: Interface binding energy of the complex
        """

        self.iam.apply(pose)
        return self.iam.get_interface_dG()


    def __ddg(self, new_pose: Pose):
        """Calculate the DDG for a new pose and evaluate it.

        Uses the cached score for the current pose and only recomputes the
        full-interface score for the proposed mutant. This avoids redundant
        old-score evaluation and speeds up each step.

        :param new_pose: New pose generated from the mutation step
        """

        new_score = self.__interface_dg_score(new_pose)
        old_score = self.current_score
        self.__algorithm(new_pose, new_score, old_score)


    def run(self):
        """Run the affinity optimization algorithm.

        The run function randomly mutates residues in the input structure,
        optionally relaxes the mutated structure, and evaluates the resulting
        pose using the selected scoring function. The algorithm will terminate
        early if the score does not improve for the specified number of steps.
        """

        if self.scoring_function == 'Distance' and not self.distances:
            raise AffinityOptimizerError('Must insert interactions using `.insert_interactions()` for the ' \
            'Distance scoring function')

        for i in range(self.number_steps):
            mutant_pose, res = random_mutation(self.current_pose, self.pos_to_mutate, self.quiet)

            if self.relax and (i % self.relax_every == 0):
                mutant_pose = fast_minimize_structure(mutant_pose, res)


            if not self.quiet:
                # prints scores
                print(f'current score: {self.current_score}')
                print(f'Best score: {self.best_score}')
                print()

            self.func(mutant_pose)

            if not self.quiet:
                progress = (i + 1) / self.number_steps
                filled = int(progress * 50)
                bar = "#" * filled + "-" * (50 - filled)
                print(f"[{bar}] {progress:.0%}")
                print()

            if self.max_no_improve == self.no_improve_steps:
                print('EARLY STOP')
                print(f'Score did not improve in {self.no_improve_steps} steps')
                print(f'Terminated algorithm at {i} step')
                break


    def view_interactions(self):
        """Display the amino acid interactions and their distances.
        Prints the amino acids, chains, positions, and calculated distance
        Prints the amino acids, chains, positions, and calculated distance
        for each interaction in the current structure.
        """
        # Get the current structure and first model
        # Get the current structure and first model
        structure = self.__get_structure()
        model = structure[0]
         # Loop through each interaction and display the residues and distance
         # Loop through each interaction and display the residues and distance
        for num, (chain, pos) in enumerate(zip(self.chains, self.positions)):
            aa1 = model[chain[0]][pos[0]].get_resname()
            aa2 = model[chain[1]][pos[1]].get_resname()
            print(f'Interaction {num + 1}')
            print(f'AA: {aa1}, chain: {chain[0]}, pos: {pos[0]}')
            print(f'AA: {aa2}, chain: {chain[1]}, pos: {pos[1]}')
            print(f'Distance: {self.distances[num]}')


    def __get_structure(self, new_pose=None):
        """Convert a PyRosetta pose into a Biopython structure.
        Uses a temporary PDB file to convert the PyRosetta pose into a
        Uses a temporary PDB file to convert the PyRosetta pose into a
        Biopython Structure object.
        """
        # if no pose ios given use current pose
        # if no pose ios given use current pose
        if new_pose is None:
            new_pose = self.current_pose
        # Create a Biopython PDB parser
        # Create a Biopython PDB parser
        parser = PDBParser(QUIET=True)

        # Create a temporary PDB file to store the PyRosetta pose
        with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
            # Write the pose to the temporary PDB file
            new_pose.dump_pdb(tmp.name)

            # Parse the temporary PDB file into a Biopython Structure object
            structure = parser.get_structure(
                "pose",
                tmp.name
            )
        return structure 


    def report_score(self):
        """Return a consistent report of optimization results.

        Returns a dict with at least the keys `best_score`, and
        `best_pose_relaxed_energy`. For DDG, the score is the interface ΔG of the
        relaxed best pose. For Distance, the score is the distance-based score of
        the relaxed best pose.
        """
        scorefxn = get_fa_scorefxn()
        relaxed_score = None
        best_pose_relaxed_energy = None
        best_distances = None
        original_score = None

        try:
            relaxed_pose = relax_structure(self.best_pose.clone())
            best_pose_relaxed_energy = scorefxn(relaxed_pose)

            if self.scoring_function == 'DDG':
                relaxed_score = self.__interface_dg_score(relaxed_pose)
                original_score = self.__interface_dg_score(self.original_pose)
            else:
                best_distances = self.find_distances(relaxed_pose)
                original_distances = self.find_distances(self.original_pose)
                relaxed_score = self.__distance_score(best_distances)
                original_score = self.__distance_score(original_distances)

        except Exception as exc:
            relaxed_pose = None
            best_pose_relaxed_energy = None
            relaxed_score = getattr(self, "best_score", None)
            print(f"report_score error: {exc}")

        report = {
            "best_score": relaxed_score,
            "best_pose_energy": best_pose_relaxed_energy,
            "scoring_function": self.scoring_function,
            "best_pose": relaxed_pose,
            "original_score": original_score

        }

        if self.scoring_function == 'Distance':
            report["best_distances"] = best_distances
            report['original_distances'] = original_distances

        return report


# endregion


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


def perform_mutation(pose, pos, amino, quiet: bool = False):
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

    if not quiet:
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


def random_mutation(pose, list_aa_pos: list, quiet: bool = False):
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
    new_pose = perform_mutation(pose, pos, aa, quiet)

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


# region RELAXING

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
        target_selector, 20.0, include_focus_in_subset=True
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

    # 4. Perform Relaxation
    relax.apply(testPose)



    return testPose

# endregion


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
