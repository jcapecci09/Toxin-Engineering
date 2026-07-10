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

    #Set Reference Pose
    RMSD_calc = pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric(pose)
    
    
    # Select Mutate Position
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))

    # Select Neighbor Position
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))
    
    # Select No Design Area
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))

    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    #True indicates here that we are flipping the selection.  So that we are turning off everything but the CDR and its neighbors.
    prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
    tf.push_back(prevent_subset_repacking)

    # Disable design
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

    # Enable design
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))
    
    # Create Packer
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    packer.task_factory(tf)
    print(tf.create_task_and_apply_taskoperations(pose))


    
    #Perform The Move
    if not os.getenv("DEBUG"):
      packer.apply(pose)


def relax_structure(pose_to_relax, output_name):

    # Create new pose and assign pose to relax to it
    testPose = Pose()
    testPose.assign(pose_to_relax)
    print(testPose)

    # Set up relax parameter
    scorefxn = get_fa_scorefxn()
    relax = rosetta.protocols.relax.FastRelax()
    relax.set_scorefxn(scorefxn)
    relax.constrain_relax_to_start_coords(True)
    print(relax)
    relax.apply(testPose)

    #rename to your desired relaxed structure name
    testPose.dump_pdb(output_name)

def perform_mutation(pdb, pos, amino):
    relaxPose = pose_from_pdb(pdb)
    # Clone it
    original = relaxPose.clone()
    scorefxn = get_score_function()


    # #Input the residue number that you wish to mutate and the 1-letter code 
    # #If you are substituting multiple, you can just have them listed separately as exampled below
    pack(relaxPose, pos, amino, scorefxn)
    # #pack(relaxPose, 81, 'D', scorefxn)
    # print("\nNew Energy:", scorefxn(relaxPose),"\n")

    original_aa = str(original.residue(pos))
    mutated_aa = relaxPose.residue(pos)

    # #SAVE THE NEW PDB FILE HERE:
    path = f'pdb_files/7K18_{pos}{amino}.pdb'
    relaxPose.dump_pdb(path)


    # #Set relaxPose back to original 
    # relaxPose = original.clone()

    print()
    print('-' * 50)
    print(f'Mutated structure succesfully saved as {path}')
    print(f'Successfully mutated {original.residue(pos).name()} at position {pos} to {relaxPose.residue(pos).name()}')
    print(f'Orginal Energy {scorefxn(original)}; New energy: {scorefxn(relaxPose)}')
    print('-' * 50)
    print('\n')


def visualize(pdb: str, output: str):
    cmd.delete("all") # Resets cell


    cmd.load(pdb)
    cmd.hide("everything")
    cmd.select('toxin', 'chain B')
    cmd.select('NaV1.5', 'chain A')
    cmd.select("binding_site", "chain A and resi 1610-1615")
    cmd.select('Histidines', 'chain B and resi 15+43')

    cmd.color('purple', 'toxin')
    cmd.color('cyan', 'NaV1.5')
    util.cbag('binding_site')
    util.cbay('Histidines')
    cmd.show('sticks', 'Histidines')
    cmd.show('sticks', 'binding_site')
    cmd.show("cartoon")

    # Center the protein
    cmd.orient()
    cmd.bg_color("white")

    cmd.png(f'Data/{output}')
    cmd.save(f"Sessions/{output}.pse")
    display(Image(f'Data/{output}.png'))