# encoding: utf-8
# env: protein_binder_design


import os
from Bio import SeqIO

from pyrosetta import init, pose_from_pdb, get_score_function
from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    InitializeFromCommandline,
    RestrictToRepacking,
    PreventRepacking,
    OperateOnResidueSubset)
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector,
    NotResidueSelector)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pose import PDBInfo



def mutate_chain_with_sequence_pyrosetta(
    pdb_filepath: str,
    chain_id: str,
    target_aa_sequence: str,
    output_pdb_filepath: str
):
    """
    Mutates a specified chain in a PDB file to a target amino acid sequence using PyRosetta.
    Removes remarks and score information from the output PDB.

    Args:
        pdb_filepath (str): Path to the input PDB file.
        chain_id (str): The ID of the chain to mutate (e.g., 'A').
        target_aa_sequence (str): The target amino acid sequence (one-letter codes).
        output_pdb_filepath (str): Path to save the mutated PDB file.
    """
    if not os.path.exists(pdb_filepath):
        print(f"Error: Input PDB file not found at {pdb_filepath}")
        return

    # PyRosetta initialization should ideally happen once for the entire script
    # and not inside a loop if calling this function multiple times.
    # For a standalone call, it's fine here.
    # For a wrapper, consider moving init() outside if it's called many times.

    print(f"Loading PDB: {pdb_filepath}")
    try:
        pose = pose_from_pdb(pdb_filepath)
    except RuntimeError as e:
        print(f"Error loading PDB into PyRosetta pose: {e}")
        return

    pdb_info = pose.pdb_info()
    chain_rosetta_indices = []
    original_chain_sequence = []

    for i in range(1, pose.total_residue() + 1):
        if pdb_info.chain(i) == chain_id and pose.residue(i).is_protein():
            chain_rosetta_indices.append(i)
            original_chain_sequence.append(pose.residue(i).name1())

    if not chain_rosetta_indices:
        print(f"Error: Chain '{chain_id}' not found or contains no protein residues in {pdb_filepath}.")
        return

    original_chain_sequence_str = "".join(original_chain_sequence)
    print(f"Original sequence of chain {chain_id} (length {len(original_chain_sequence_str)}): {original_chain_sequence_str}")

    if len(chain_rosetta_indices) != len(target_aa_sequence):
        print(f"Error: Length of target sequence ({len(target_aa_sequence)}) "
              f"does not match length of chain {chain_id} protein residues ({len(chain_rosetta_indices)}).")
        return

    print(f"Target sequence for chain {chain_id} (length {len(target_aa_sequence)}): {target_aa_sequence}")

    pose_to_mutate = pose.clone()
    mutations_applied_count = 0
    print("Applying mutations...")
    for i, rosetta_idx in enumerate(chain_rosetta_indices):
        original_aa = pose_to_mutate.residue(rosetta_idx).name1()
        target_aa = target_aa_sequence[i].upper()

        if original_aa != target_aa:
            try:
                mutate_residue(pose_to_mutate, rosetta_idx, target_aa, pack_radius=0.0)
                mutations_applied_count += 1
            except RuntimeError as e:
                print(f"  Error mutating chain {chain_id} Ros_idx {rosetta_idx} from {original_aa} to {target_aa}: {e}")
                print(f"  This might be due to an invalid target AA code ('{target_aa}') or other Rosetta issue. Skipping.")
                continue

    if mutations_applied_count > 0:
        print(f"Applied {mutations_applied_count} mutations to chain {chain_id}.")
    elif original_chain_sequence_str.upper() == target_aa_sequence.upper():
        print(f"Target sequence is identical to the original sequence of chain {chain_id}. No mutations applied.")
    else:
        print(f"No mutations were effectively applied, but sequences may differ. Check for errors above.")

    print(f"Optimizing structure around mutated chain {chain_id}...")
    scorefxn = get_score_function(True)

    chain_selector = ChainSelector(chain_id)

    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())
    tf.push_back(RestrictToRepacking())
    # prevent_repacking_operation = PreventRepacking()
    # not_target_chain_selector = NotResidueSelector(chain_selector)
    # tf.push_back(OperateOnResidueSubset(prevent_repacking_operation, not_target_chain_selector))

    print("Performing localized sidechain repacking...")
    packer = PackRotamersMover(scorefxn, tf)
    try:
        packer.apply(pose_to_mutate)
        print(f"Sidechain repacking completed for chain {chain_id}.")
    except RuntimeError as e:
        print(f"Error during sidechain repacking for chain {chain_id}: {e}")
        print("Proceeding, but sidechains might not be optimal.")


    try:
        pose_to_mutate.dump_pdb(output_pdb_filepath)
        print(f"Successfully saved mutated PDB to: {output_pdb_filepath}")
    except RuntimeError as e:
        print(f"Error saving mutated PDB: {e}")




def process_fasta_for_mutation(
    fasta_filepath: str,
    input_pdb_filepath: str,
    chain_id: str,
    output_directory: str
):
    """
    Reads a FASTA file, and for each record, mutates the specified chain in the
    input PDB file to the sequence from the FASTA record.
    Output PDBs are named after the FASTA record IDs.

    Args:
        fasta_filepath (str): Path to the input FASTA file containing target sequences.
        input_pdb_filepath (str): Path to the base PDB file to be mutated.
        chain_id (str): The ID of the chain in the PDB to mutate (e.g., 'A').
        output_directory (str): Directory to save the mutated PDB files.
    """
    if not os.path.exists(fasta_filepath):
        print(f"Error: FASTA file not found at {fasta_filepath}")
        return

    if not os.path.exists(input_pdb_filepath):
        print(f"Error: Input PDB file not found at {input_pdb_filepath}")
        return

    os.makedirs(output_directory, exist_ok=True)
    print(f"Output directory set to: {output_directory}")

    print(f"\nInitializing PyRosetta for the session...")
    # init_options = (
    #     "-ignore_unrecognized_res true "
    #     "-load_PDB_components false "
    #     "-mute core.conformation.Conformation core.scoring core.pack "
    #     "-ex1 -ex2 "
    #     "-no_optH false " # Ensure hydrogens are present/optimized
    #     "-out:file:no_score_header true " # Suppress score line at the top
    #     "-out:file:no_score_line true " # Suppress score line at the bottom
    #     "-out:file:no_timestamp true " # Suppress timestamp remarks
    #     "-out:file:no_comments true " # Suppress general comments/remarks
    # )

    init_options = (
        "-ignore_unrecognized_res true "
        "-load_PDB_components false "
        "-mute core.conformation.Conformation core.scoring core.pack "
        "-ex1 -ex2"
    )

    try:
        # Initialize PyRosetta once for the entire session
        init(options=init_options, silent=False)
    except Exception as e:
        print(f"PyRosetta initialization failed: {e}")
        print("This can happen if PyRosetta is already initialized with different options.")
        print("Please restart your Python session if you encounter this frequently.")
        return

    print(f"\nProcessing FASTA file: {fasta_filepath}")
    processed_count = 0
    try:
        for record in SeqIO.parse(fasta_filepath, "fasta"):
            fasta_id = record.id.replace(os.sep, "_").replace(" ", "_") # Sanitize ID for filename
            target_sequence = str(record.seq).upper() # Ensure uppercase

            output_pdb_filename = f"{fasta_id}.pdb"
            output_pdb_filepath = os.path.join(output_directory, output_pdb_filename)

            print(f"\n--- Processing record: {fasta_id} ---")
            print(f"  Target sequence: {target_sequence}")
            print(f"  Saving to: {output_pdb_filepath}")

            mutate_chain_with_sequence_pyrosetta(
                pdb_filepath=input_pdb_filepath,
                chain_id=chain_id,
                target_aa_sequence=target_sequence,
                output_pdb_filepath=output_pdb_filepath
            )
            processed_count += 1
            print(f"--- Finished processing record: {fasta_id} ---")

    except Exception as e:
        print(f"Error reading FASTA file or during processing: {e}")
        print("Please ensure Biopython is installed (`pip install biopython`) and the FASTA file is correctly formatted.")
        return

    print(f"\n--- Processing complete ---")
    print(f"Processed {processed_count} sequences from {fasta_filepath}.")



if __name__ == "__main__":

    fasta_file = "/media/hao/Data/Project_kinase_substrate_design/2025-05-12_ALK_sub/previous_sub.fasta"
    input_pdb = "/media/hao/Data/Project_kinase_substrate_design/2025-05-12_ALK_sub/met_backbone.pdb"
    target_chain = "A" # The chain in your PDB to mutate
    output_dir = "/media/hao/Data/Project_kinase_substrate_design/2025-05-12_ALK_sub/met_subs"

    # mutate_chain_with_sequence_pyrosetta(
    #     pdb_filepath=input_pdb,
    #     chain_id=target_chain,
    #     target_aa_sequence="GDDDEYVDDES",  # Example sequence, replace with actual
    #     output_pdb_filepath=os.path.join(output_dir, "mutated_chain.pdb")  # Example output path
    # )
    

    process_fasta_for_mutation(fasta_file, input_pdb, target_chain, output_dir)