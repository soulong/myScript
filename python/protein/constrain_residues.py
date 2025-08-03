# -*- coding: utf-8 -*-
# env: proteinmpnn_binder_design

import argparse
from Bio.PDB import PDBParser, PDBIO
import os
import sys

# Mapping from single-letter to three-letter amino acid codes
SINGLE_TO_THREE_LETTER_AA = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

def constrain_residues(input_pdb_file, constraints, output_pdb_file):
    """
    Modifies a PDB file to constrain specific residues based on absolute positions
    or relative middle positions.

    Args:
        input_pdb_file (str): Path to the input PDB file.
        constraints (list): A list of dictionaries, where each dictionary represents
                            a constraint with keys: 'chain', 'res_num', 'aa', 'type'.
                            'type' can be 'absolute' or 'middle'.
        output_pdb_file (str): Path to save the modified PDB file.
    Returns:
        list: A list of (chain_id, residue_id) tuples for all residues that were constrained.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("peptide", input_pdb_file)
    except Exception as e:
        print(f"Error parsing {input_pdb_file}: {e}")
        return []

    io = PDBIO()
    io.set_structure(structure)

    # Store indices of changed residues for remarks
    changed_residue_indices = [] # Stores (chain_id, residue_id) tuples

    # Group constraints by chain for efficient processing
    chain_specific_constraints = {}
    for constraint in constraints:
        chain_id = constraint['chain']
        if chain_id not in chain_specific_constraints:
            chain_specific_constraints[chain_id] = []
        chain_specific_constraints[chain_id].append(constraint)


    for model in structure:
        for chain in model:
            # Process only chains that have constraints defined
            if chain.id in chain_specific_constraints:
                residues = list(chain)
                chain_length = len(residues)

                for constraint in chain_specific_constraints[chain.id]:
                    con_res_num = constraint['res_num']
                    con_aa_3letter = constraint['aa']
                    con_type = constraint['type']

                    if con_type == 'absolute':
                        found = False
                        for residue in chain:
                            if residue.id[1] == con_res_num:
                                residue.resname = con_aa_3letter
                                changed_residue_indices.append((chain.id, residue.id[1]))
                                found = True
                                break
                        if not found:
                            print(f"Warning: Absolute residue {con_res_num} not found in chain {chain.id} of {input_pdb_file}. Skipping.")
                    elif con_type == 'middle':
                        if chain_length == 0:
                            print(f"Warning: Chain {chain.id} in {input_pdb_file} has zero length. Cannot apply middle constraint.")
                            continue

                        middle_index = chain_length // 2
                        target_index = middle_index + con_res_num # con_res_num is now the relative offset

                        if 0 <= target_index < chain_length:
                            residue_to_constrain = residues[target_index]
                            residue_to_constrain.resname = con_aa_3letter
                            changed_residue_indices.append((chain.id, residue_to_constrain.id[1]))
                        else:
                            print(f"Warning: Relative middle position {con_res_num} out of bounds for chain {chain.id} in {input_pdb_file} (length: {chain_length}). Skipping.")

    try:
        io.save(output_pdb_file)
    except Exception as e:
        print(f"Error saving {output_pdb_file}: {e}")
        return []

    return changed_residue_indices


def clean_pdb_rfdiffusion(pdb_path):
    """
    Reads a PDB file, removes TER and END records, and strips the Element symbol (cols 77-78)
    and Charge (cols 79-80) columns from ATOM and HETATM records.
    The file is overwritten with the cleaned content.
    This cleaning is typical for preparing PDBs for tools like RFdiffusion.

    Args:
        pdb_path (str): Path to the PDB file to clean.
    """
    temp_lines = []
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                # Remove TER and END records
                if line.startswith("TER") or line.strip() == "END":
                    continue

                if line.startswith(("ATOM  ", "HETATM")):
                    # PDB ATOM/HETATM record column definitions:
                    # 1-6  Record name "ATOM  " or "HETATM"
                    # ...
                    # 73-76  LString(4)  segID Segment identifier
                    # 77-78  LString(2)  element Element symbol (right-justified)
                    # 79-80  LString(2)  charge   Charge on the atom (right-justified)

                    # Keep everything up to column 76 (index 75)
                    cleaned_line = line[0:76] # Slice up to, but not including, index 76
                    temp_lines.append(cleaned_line + '\n') # Add newline back
                else:
                    temp_lines.append(line) # Keep other lines as is
        
        # # Ensure a single END record at the very end of the cleaned file
        # # This handles cases where the original END was removed or multiple ENDs existed.
        # if temp_lines and not temp_lines[-1].strip() == 'END':
        #     temp_lines.append('END\n')
        # elif not temp_lines: # If file became empty after cleaning, just add END
        #     temp_lines.append('END\n')

        with open(pdb_path, 'w') as f:
            f.writelines(temp_lines)
        print(f"Cleaned PDB for RFdiffusion in {pdb_path} (removed TER/END records and element/charge columns).")

    except Exception as e:
        print(f"Error cleaning PDB for RFdiffusion in {pdb_path}: {e}")



def add_remarks_to_pdb(pdb_path, changed_residue_indices):
    """
    Adds REMARK lines to a PDB file to indicate constrained residues.
    The remarks are inserted before the END record.

    Args:
        pdb_path (str): Path to the PDB file to modify.
        changed_residue_indices (list): A list of (chain_id, residue_id) tuples for constrained residues.
    """
    if not changed_residue_indices:
        return # No remarks to add

    remarks_to_add = []
    # Sort for consistent output
    sorted_indices = sorted(changed_residue_indices, key=lambda x: (x[0], x[1]))
    for chain_id, res_id in sorted_indices:
        # PDB format for REMARK is fixed at 79 chars, usually REMARK + 3 spaces + content
        remark = f"REMARK PDBinfo-LABEL:{res_id: >5} FIXED"
        remarks_to_add.append(remark)

    try:
        remarks_str = '\n'.join(remarks_to_add)
        with open(pdb_path, 'a') as f:
            f.write('\n')
            f.write(remarks_str)
        print(f"Added remarks to {pdb_path}.")
        
    except Exception as e:
        print(f"Error adding remarks to {pdb_path}: {e}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Constrain residues in PDB files based on absolute or relative middle positions, using a unified input format and single-letter AA codes, and add remarks to output PDBs.")
    parser.add_argument("input_pdb", help="Path to the input PDB file or a directory containing PDB files.")
    parser.add_argument("--constraints", help='Comma-separated list of constraints in the format "chain:resnum_or_relpos:AA". '
                                             'For relative middle positions, prefix resnum with "M" (e.g., "A:M0:A"). '
                                             'For absolute positions, use just the number (e.g., "B:10:P"). '
                                             'AA should be a single-letter amino acid code.',
                        default="")
    parser.add_argument("--output_dir", help="Directory to save the modified PDB files. If not specified, overwrites the input file.", default=None)
    parser.add_argument("--suffix", help="Suffix to add to the output file names.", default="")

    args = parser.parse_args()

    all_constraints = []
    if args.constraints:
        for constraint_str in args.constraints.split(','):
            parts = constraint_str.split(':')
            if len(parts) == 3:
                chain_id, pos_str, single_letter_aa = parts
                
                if not single_letter_aa.upper() in SINGLE_TO_THREE_LETTER_AA:
                    print(f"Error: Invalid single-letter amino acid code '{single_letter_aa}' in constraint '{constraint_str}'.")
                    sys.exit(1)
                three_letter_aa = SINGLE_TO_THREE_LETTER_AA[single_letter_aa.upper()]

                try:
                    if pos_str.startswith('M'):
                        # Middle constraint: 'M' followed by relative position
                        relative_position = int(pos_str[1:])
                        all_constraints.append({
                            'chain': chain_id.strip(),
                            'res_num': relative_position, # this is the relative offset
                            'aa': three_letter_aa,
                            'type': 'middle'
                        })
                    else:
                        # Absolute constraint: just the residue number
                        absolute_residue_number = int(pos_str)
                        all_constraints.append({
                            'chain': chain_id.strip(),
                            'res_num': absolute_residue_number,
                            'aa': three_letter_aa,
                            'type': 'absolute'
                        })
                except ValueError:
                    print(f"Error: Invalid position format '{pos_str}' in constraint '{constraint_str}'. Must be an integer or 'M' followed by an integer.")
                    sys.exit(1)
            else:
                print(f"Error: Invalid constraint format '{constraint_str}'. Expected 'chain:resnum_or_relpos:AA'.")
                sys.exit(1)


    if os.path.isfile(args.input_pdb):
        input_files = [args.input_pdb]
    elif os.path.isdir(args.input_pdb):
        input_files = [os.path.join(args.input_pdb, f) for f in os.listdir(args.input_pdb) if f.endswith(".pdb")]
    else:
        print(f"Error: Input path '{args.input_pdb}' is not a file or directory.")
        sys.exit(1)

    for input_file in input_files:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        if args.output_dir:
            output_file = os.path.join(args.output_dir, f"{base_name}{args.suffix}.pdb")
            os.makedirs(args.output_dir, exist_ok=True)
        else:
            output_file = input_file # Overwrite original if no output_dir specified

        print(f"Processing: {input_file} -> {output_file}")
        
        # Step 1: Apply constraints and get list of changed residues
        changed_res_info = constrain_residues(input_file, all_constraints, output_file)
        # Step 2: Optionally clean PDB for RFdiffusion
        clean_pdb_rfdiffusion(output_file)
        # Step 3: Add remarks to the saved PDB file
        add_remarks_to_pdb(output_file, changed_res_info)


    print("Finished processing PDB files.")