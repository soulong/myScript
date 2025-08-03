#!/usr/bin/env python3
"""
Extract sequence IDs from FASTA file and write to a text file.
Output file name is based on input FASTA file name unless specified.
"""

import argparse
from pathlib import Path

def extract_fasta_ids(fasta_file, output_file=None):
    """
    Extract sequence IDs from a FASTA file and write to text file.
    
    Args:
        fasta_file (str): Path to FASTA file
        output_file (str): Output text file path (optional)
    """
    fasta_path = Path(fasta_file)
    
    if not fasta_path.exists():
        print(f"Error: FASTA file {fasta_file} does not exist.")
        return
    
    # Generate output file name if not specified
    if output_file is None:
        output_file = fasta_path.with_suffix('.txt')
    else:
        output_file = Path(output_file)
    
    sequence_ids = []
    
    try:
        with open(fasta_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line.startswith('>'):
                    # Extract sequence ID (everything after '>')
                    seq_id = line[1:]  # Remove the '>' character
                    
                    # Handle cases where there might be description after the ID
                    # Take only the first part before any whitespace
                    seq_id = seq_id.split()[0] if seq_id.split() else seq_id
                    
                    sequence_ids.append(seq_id)
        
        if not sequence_ids:
            print(f"No sequence IDs found in {fasta_file}")
            return
        
        # Write sequence IDs to output file
        with open(output_file, 'w') as f:
            for seq_id in sequence_ids:
                f.write(seq_id + '\n')
        
        print(f"Extracted {len(sequence_ids)} sequence IDs from {fasta_path.name}")
        print(f"Output written to: {output_file}")
        
        # Show first few IDs as preview
        preview_count = min(5, len(sequence_ids))
        print(f"\nFirst {preview_count} sequence IDs:")
        for i, seq_id in enumerate(sequence_ids[:preview_count]):
            print(f"  {i+1}. {seq_id}")
        
        if len(sequence_ids) > preview_count:
            print(f"  ... and {len(sequence_ids) - preview_count} more")
    
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Extract sequence IDs from FASTA file')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('-o', '--output', help='Output text file (default: same name as FASTA with .txt extension)')
    
    args = parser.parse_args()
    
    extract_fasta_ids(args.fasta_file, args.output)

if __name__ == "__main__":
    main()