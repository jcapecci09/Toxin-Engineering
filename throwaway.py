#Python
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.teaching import *
from pyrosetta.toolbox import *

#Core Includes
#from rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *

#Protocol Includes
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax

# relevant for creating pymol session and images
from pymol import cmd
from IPython.display import Image, display
from pymol import util
import tempfile

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
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)

    # Select Neighbor Position
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_distance(8.0) # set neighborhood distance to 8 angstroms
    nbr_selector.set_focus_selector(mut_posi) # Build neighborhood around residue
    nbr_selector.set_include_focus_in_subset(True) # include mutant residue in neighborhood
    
    # Select No Design Area
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)

    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These are TaskOperations which are Instructions for packing

    # Intialize rules stated when you init()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())

    # Includes the current side-chain conformation (rotamer) as one of the candidates
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent()) 

    # Prevents residues involved in disulfide bonds (cysteine–cysteine bonds) from being repacked.
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    #True indicates here that we are flipping the selection.  So that we are turning off everything but the CDR and its neighbors.
    prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
    tf.push_back(prevent_subset_repacking)

    # Disable design
    # only change identity of mutant AA
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), not_design))

    # Enable design
    # allows amino acid idnetity to change
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))
    
    # Create Packer
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    packer.task_factory(tf)

    # Perform The Move
    if not os.getenv("DEBUG"):
      packer.apply(pose)


def relax_structure(pose_to_relax, output_name):
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
    relax = rosetta.protocols.relax.FastRelax()

    # Use the full-atom score function during relaxation
    relax.set_scorefxn(scorefxn)

    # Keep the relaxed structure close to the starting coordinates
    relax.constrain_relax_to_start_coords(True)

    # Perform energy minimization and side-chain optimization
    relax.apply(testPose)

    # Print energy after relaxed
    print(f"Relaxed energy: {scorefxn(testPose)}")

    # rename to your desired relaxed structure name

    if output_name is not None:
        testPose.dump_pdb(output_name)
    else:
        return testPose


def perform_mutation(pose_to_mutate, pos, amino):
    """Wrapper function to perform pack without specifying scoring function

    :param relaxed_pdb: pdb to perform mutation on
    :param pos: position to perform mutation
    :param amino: amino acid to muatate in
    :return: pose of muatated pdb
    """

    # Clone it
    original = pose_to_mutate.clone()

    # Create the default Rosetta score function
    scorefxn = get_score_function()

    # Perform packing
    pack(pose_to_mutate, pos, amino, scorefxn)


    # print useful information
    print()
    print('-' * 50)
    print(f'TASK COMPLETE')
    print(f'Successfully mutated {original.residue(pos).name()} at position {pos} to {pose_to_mutate.residue(pos).name()}')
    print(f'Orginal Energy {scorefxn(original)}; New energy: {scorefxn(pose_to_mutate)}')
    print('-' * 50)
    print('\n')

    return pose_to_mutate

def visualize(pdb: str, output: str):
    """Allows protein complex to be viewed in PyMol. Outputs both
    a pymol session in the Sessions directory and a quick visualization. 

    :param pdb: pdb to view
    :param output: output path to place PyMol session
    """

    # Reset cell to ensure no other proteins are in the current session
    cmd.delete("all")

    # perform basic tasks
    cmd.load(pdb) # load pdb file
    cmd.hide("everything") # hide everything

    # select the toxin and color it purple
    cmd.select('toxin', 'chain B') 
    cmd.color('purple', 'toxin')

    # Select the sodium channel and color it cyan
    cmd.select('NaV1.5', 'chain A')
    cmd.color('cyan', 'NaV1.5')

    # Select the binding site on the sodium channel and color it
    cmd.select("binding_site", "chain A and resi 1610-1615")
    util.cbag('binding_site')

    # Select the toxins histidines and color them
    cmd.select('Histidines', 'chain B and resi 15+43')
    util.cbay('Histidines')

    # Show sticks for the relevant amino acids
    cmd.show('sticks', 'Histidines')
    cmd.show('sticks', 'binding_site')

    # Show cartoon for whole complex
    cmd.show("cartoon")

    # Center the protein around interface
    cmd.set_view((\
     0.843022108,    0.219623774,   -0.490965903,\
    -0.473968744,    0.734833062,   -0.485122830,\
     0.254241019,    0.641686916,    0.723597884,\
     0.002971664,   -0.004409352,  -47.775711060,\
   102.082580566,  132.507263184,  167.430847168,\
  -1443.574340820, 1536.576416016,  -19.999998093 ))

    # Set background to white
    cmd.bg_color("white")

    # Save png and session
    cmd.png(f'Data/{output}')
    cmd.save(f"Sessions/{output}.pse")

    # Display png for viewing in notebook
    display(Image(f'Data/{output}.png'))


def run_analysis(relaxed_pdb, pos, amino_acid): 
    """Wrapper to perform entire mutation analysis. Mutate a pdb file, 
    relax it, and then visualize in PyMol

    :param relaxed_pdb: pdgb file that has already been relaxed by rosetta
    :param pos: postion to mutate
    :param amino_acid: amino acid to muatate in
    """

    # Grab the output path
    output_path = f'{relaxed_pdb.split('/')[1].split('.')[0]}_{pos}{amino_acid}'

    # perform mutation
    pose = perform_mutation(relaxed_pdb, pos, amino_acid)

    # Relax the structure
    relax_structure(pose, f'pdb_files/{output_path}.pdb')

    # View in PyMol
    visualize(f'pdb_files/{output_path}.pdb', f'{output_path}')



def insert_mutation(pdb, mutation, seq_mutated):

    pose = pose_from_pdb(pdb)

    new_pose = Pose()
    new_pose.assign(pose)

    pose_seq = pose.sequence()

    start = pose_seq.find(seq_mutated)

    for pos, aa in zip(range(start, len(mutation) + start), mutation):
        new_pose = relax_structure(perform_mutation(new_pose, pos + 1, aa), None)

    return new_pose





# pose = insert_mutation('coding/Toxin-Engineering/pdb_files/7K18_relax.pdb', 'KVCGCDLAQGGFF', 'GTVLSDIIQKYFF')

# pose.dump_pdb('coding/Toxin-Engineering/pdb_files/7K18_mutant.pdb')wwwwwwa d d

def seperate_chain(pose):
    begin = pose.conformation().chain_begin(1)
    end = pose.conformation().chain_end(1)
    begin2 = pose.conformation().chain_begin(2)
    end2 = pose.conformation().chain_end(2)

    chain1 = pose.clone()
    chain2 = pose.clone()

    delete_region(chain1, begin, end)
    delete_region(chain2, begin2, end2)

    return chain1, chain2