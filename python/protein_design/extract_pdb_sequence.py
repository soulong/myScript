#!/usr/bin/env python3
"""
Extract chain A sequences from all PDB files in a directory and write to FASTA format.
Uses BioPython for robust PDB parsing and sequence extraction.
"""

import os
import sys
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import argparse

# Standard amino acid three-letter to one-letter code mapping
AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def extract_chain_sequence(pdb_file, chain_id='A'):
    """
    Extract amino acid sequence from a specific chain in a PDB file.
    
    Args:
        pdb_file (str): Path to PDB file
        chain_id (str): Chain identifier (default: 'A')
    
    Returns:
        str: Amino acid sequence in single-letter code, or None if chain not found
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('protein', pdb_file)
        
        # Get the first model (most PDB files have only one model)
        model = structure[0]
        
        # Check if chain exists
        if chain_id not in model:
            return None
        
        chain = model[chain_id]
        sequence = []
        
        # Extract residues from chain
        for residue in chain:
            # Skip non-amino acid residues (water, ligands, etc.)
            if residue.id[0] == ' ':  # Standard amino acid residues have ' ' as hetfield
                res_name = residue.get_resname()
                if res_name in AA_MAP:
                    sequence.append(AA_MAP[res_name])
                else:
                    # Handle non-standard amino acids by replacing with 'X'
                    sequence.append('X')
        
        return ''.join(sequence)
    
    except Exception as e:
        print(f"Error processing {pdb_file}: {e}")
        return None



def process_pdb_input(input_path, chain_id='A', output_file='chain_sequences.fasta'):
    """
    Process PDB input - can be a single file or directory containing PDB files.
    
    Args:
        input_path (str): Path to PDB file or directory containing PDB files
        chain_id (str): Chain identifier to extract (default: 'A')
        output_file (str): Output FASTA file name
    """
    input_path = Path(input_path)
    
    if not input_path.exists():
        print(f"Error: {input_path} does not exist.")
        return
    
    pdb_files = []
    
    if input_path.is_file():
        # Single PDB file
        if input_path.suffix.lower() == '.pdb':
            pdb_files = [input_path]
        else:
            print(f"Error: {input_path} is not a PDB file (.pdb extension required)")
            return
    elif input_path.is_dir():
        # Directory containing PDB files
        pdb_files = list(input_path.glob('*.pdb')) + list(input_path.glob('*.PDB'))
        if not pdb_files:
            print(f"No PDB files found in {input_path}")
            return
    else:
        print(f"Error: {input_path} is neither a file nor a directory.")
        return
    
    sequences = []
    processed = 0
    skipped = 0
    
    print(f"Processing {len(pdb_files)} PDB file(s)...")
    
    for pdb_file in sorted(pdb_files):
        # Get filename without extension as sequence ID
        seq_id = pdb_file.stem
        
        # Extract sequence from specified chain
        sequence = extract_chain_sequence(pdb_file, chain_id)
        
        if sequence:
            sequences.append((seq_id, sequence))
            processed += 1
            print(f"✓ {seq_id}: {len(sequence)} residues")
        else:
            skipped += 1
            print(f"✗ {seq_id}: Chain {chain_id} not found or empty")
    
    # Write sequences to FASTA file
    if sequences:
        with open(output_file, 'w') as f:
            for seq_id, sequence in sequences:
                f.write(f">{seq_id}\n")
                # Write sequence in lines of 80 characters (standard FASTA format)
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')
        
        print(f"\nResults:")
        print(f"  Processed: {processed} file(s)")
        print(f"  Skipped: {skipped} file(s)")
        print(f"  Output written to: {output_file}")
    else:
        print("No sequences were extracted.")

def main():
    parser = argparse.ArgumentParser(description='Extract chain sequences from PDB file(s)')
    parser.add_argument('input', help='PDB file or directory containing PDB files')
    parser.add_argument('-c', '--chain', default='A', help='Chain ID to extract (default: A)')
    parser.add_argument('-o', '--output', default='chain_sequences.fasta', 
                       help='Output FASTA file (default: chain_sequences.fasta)')
    
    args = parser.parse_args()
    
    process_pdb_input(args.input, args.chain, args.output)

if __name__ == "__main__":
    main()