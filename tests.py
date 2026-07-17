from pyrosetta import *
import os
from pymol import cmd
import tempfile
from pymol import util

init()

def motif_pose(pose, motif):
    pose_seq = pose.sequence()

    length = len(pose_seq)
    start_res = pose_seq.find(motif) + 1
    end_res = start_res + len(motif) - 1

    pyrosetta.rosetta.protocols.grafting.delete_region(pose, end_res + 1, pose.total_residue())
    pyrosetta.rosetta.protocols.grafting.delete_region(pose, 1, start_res - 1)



def structural_comparison(pose1,  pose2, motif):

    # Create copies to ensure original poses aren't
    # being altered
    trunc_pose1 = Pose()
    trunc_pose2 = Pose()
    trunc_pose1.assign(pose1)
    trunc_pose2.assign(pose2)

    motif_pose(trunc_pose1, motif)
    print(trunc_pose1.sequence())
    motif_pose(trunc_pose2, motif)
    print(trunc_pose2.sequence())


    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp1, tempfile.NamedTemporaryFile(suffix=".pdb") as tmp2:
    
        trunc_pose1.dump_pdb(tmp1.name)
        trunc_pose2.dump_pdb(tmp2.name)

        cmd.load(tmp1.name, "pose1")
        cmd.load(tmp2.name, "pose2")

        rmsd = cmd.align("pose2 and name N+CA+C+O", "pose1 and name N+CA+C+O")
        cmd.show("cartoon", "pose1")
        cmd.show("cartoon", "pose2")

        cmd.show("sticks", "pose1")
        cmd.show("sticks", "pose2")

        util.cbag('pose1')
        util.cbay('pose2')
        cmd.save("coding/Toxin-Engineering/Sessions/alignment_session.pse")

        print("RMSD:", rmsd)


    return rmsd

pose1 = pose_from_pdb('/home/jcape/coding/Toxin-Engineering/pdb_files/7K18_mutant.pdb')
pose2 = pose_from_pdb('/home/jcape/coding/Toxin-Engineering/ABLIM1_alpha.pdb')
structural_comparison(pose1, pose2, 'KVCGCDLAQGGFF')

