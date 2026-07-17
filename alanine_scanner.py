import numpy as np
from Bio.PDB import PDBParser
import pandas as pd
from utils import perform_mutation
from pyrosetta import get_fa_scorefxn, pose_from_pdb
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover


def contacts(pdb: str, cutoff: int) -> tuple[set[int], set[int]]:


    protein_chain_A = "A"
    protein_chain_B = "B"

    # ============================================================

    # READ PDB


    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb)
    model = structure[0]

    protein_A = model[protein_chain_A]
    protein_B = model[protein_chain_B]

    # Keep only standard amino acids
    protein_A_residues = [r for r in protein_A if r.id[0] == " "]
    protein_B_residues   = [r for r in protein_B if r.id[0] == " "]

    print(f"Protein_A residues : {len(protein_A_residues)}")
    print(f"Protein_B residues   : {len(protein_B_residues)}")

    # ============================================================

    # COMPUTE MINIMUM HEAVY-CHAIN ATOM DISTANCE BETWEEN EVERY RESIDUE PAIR


    contact_matrix = np.zeros((len(protein_A_residues),
                            len(protein_B_residues)))

    for i, A_res in enumerate(protein_A_residues):

        A_atoms = [a for a in A_res if a.element != "H"]

        for j, B_res in enumerate(protein_B_residues):

            B_atoms = [a for a in B_res if a.element != "H"]

            min_dist = np.inf

            for A_a in A_atoms:
                for B_a in B_atoms:
                    d = A_a - B_a
                    if d < min_dist:
                        min_dist = d

            contact_matrix[i, j] = min_dist


    plot_matrix = contact_matrix.copy()



    # Hide anything farther than cutoff
    plot_matrix[plot_matrix > cutoff] = np.nan


    # Make dataframe of contacts
    df = pd.DataFrame(plot_matrix)
    contacts = df.stack().dropna().reset_index()
    contacts.columns = ['Chain A', 'Chain B', 'Distance']
    contacts.loc[:, 'Chain B'] += 1
    contacts.loc[:, 'Chain A'] += 1

    return set(contacts.loc[:, 'Chain A']), set(contacts.loc[:, 'Chain B'])



def alanine_scanner(pdb: str, contacts: list[set, set]):

    
    pose = pose_from_pdb(pdb)
    end = pose.conformation().chain_end(1)
    pos_seq = pose.sequence()


    original_dg = delta_g(pose)

    data = {}
    for chain_num, set_of_contacts in enumerate(contacts):
        for pos in set_of_contacts:
            mutant_pose = perform_mutation(pose, pos, 'A')
            if chain_num == 1:
                aa_pos = pos + end
            else:
                aa_pos = pos
            data[pos] = {"ddG_binding": delta_g(mutant_pose) - original_dg, 
                         "chain_num": chain_num + 1, 'AA': pos_seq[aa_pos]}
            

    return pd.DataFrame(data).T
    
    

def delta_g(pose):

    scorefxn = get_fa_scorefxn()
    iam = InterfaceAnalyzerMover(
    1,          # interface jump
    False,      # tracer output
    scorefxn)
    iam.apply(pose)
    return iam.get_interface_dG()


contact = contacts('/home/jcape/coding/Toxin-Engineering/pdb_files/7K18_relax.pdb', 6)
print(alanine_scanner('/home/jcape/coding/Toxin-Engineering/pdb_files/7K18_relax.pdb', contact))

