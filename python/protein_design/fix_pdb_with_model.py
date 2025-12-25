#!/usr/bin/env python3
"""
conda env: proteinmpnn_binder_design
Script to complete missing residues in PDB structures using AlphaFold predictions.
Uses structural alignment to identify missing regions and fills them with AlphaFold data.
"""

import os
import sys
import argparse
from pathlib import Path
import numpy as np
from Bio import PDB
from Bio.PDB import PDBIO, Structure, Model, Chain, Residue, Atom
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

# Suppress PDB construction warnings
warnings.simplefilter('ignore', PDBConstructionWarning)

class PDBCompleter:
    def __init__(self, pdb_file, alphafold_file, output_file=None):
        """
        Initialize the PDB completer.
        
        Args:
            pdb_file (str): Path to the PDB structure file
            alphafold_file (str): Path to the AlphaFold structure file
            output_file (str): Path for output completed structure
        """
        self.pdb_file = pdb_file
        self.alphafold_file = alphafold_file
        self.output_file = output_file or self._generate_output_name()
        
        self.parser = PDB.PDBParser(QUIET=True)
        self.pdb_structure = None
        self.af_structure = None
        
    def _generate_output_name(self):
        """Generate output filename based on input PDB file."""
        base_name = Path(self.pdb_file).stem
        return f"{base_name}_completed.pdb"
    
    def load_structures(self):
        """Load both PDB and AlphaFold structures."""
        try:
            self.pdb_structure = self.parser.get_structure("pdb", self.pdb_file)
            self.af_structure = self.parser.get_structure("alphafold", self.alphafold_file)
            print(f"Loaded PDB structure: {self.pdb_file}")
            print(f"Loaded AlphaFold structure: {self.alphafold_file}")
        except Exception as e:
            print(f"Error loading structures: {e}")
            sys.exit(1)
    
    def get_residue_numbers(self, chain):
        """Get list of residue numbers in a chain."""
        return [res.id[1] for res in chain if res.id[0] == ' ']
    
    def find_missing_residues(self, pdb_chain, af_chain):
        """
        Find residues that are missing in PDB but present in AlphaFold.
        
        Args:
            pdb_chain: PDB chain object
            af_chain: AlphaFold chain object
            
        Returns:
            list: List of missing residue numbers
        """
        pdb_residues = set(self.get_residue_numbers(pdb_chain))
        af_residues = set(self.get_residue_numbers(af_chain))
        
        missing_residues = sorted(af_residues - pdb_residues)
        return missing_residues
    
    def get_common_residues(self, pdb_chain, af_chain):
        """Get residues common to both structures for alignment."""
        pdb_residues = set(self.get_residue_numbers(pdb_chain))
        af_residues = set(self.get_residue_numbers(af_chain))
        
        common_residues = sorted(pdb_residues & af_residues)
        return common_residues
    
    def align_structures(self, pdb_chain, af_chain):
        """
        Align AlphaFold structure to PDB structure using common residues.
        
        Args:
            pdb_chain: PDB chain object
            af_chain: AlphaFold chain object
            
        Returns:
            tuple: (superimposer, rmsd)
        """
        common_residues = self.get_common_residues(pdb_chain, af_chain)
        
        if len(common_residues) < 3:
            print(f"Warning: Only {len(common_residues)} common residues found. Alignment may be poor.")
            return None, float('inf')
        
        # Get CA atoms for alignment
        pdb_atoms = []
        af_atoms = []
        
        for res_num in common_residues:
            try:
                pdb_res = pdb_chain[res_num]
                af_res = af_chain[res_num]
                
                if 'CA' in pdb_res and 'CA' in af_res:
                    pdb_atoms.append(pdb_res['CA'])
                    af_atoms.append(af_res['CA'])
            except KeyError:
                continue
        
        if len(pdb_atoms) < 3:
            print("Error: Not enough CA atoms for alignment")
            return None, float('inf')
        
        # Perform superimposition
        superimposer = Superimposer()
        superimposer.set_atoms(pdb_atoms, af_atoms)
        
        # Apply transformation to entire AlphaFold structure
        af_atoms_all = []
        for model in self.af_structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        af_atoms_all.append(atom)
        
        superimposer.apply(af_atoms_all)
        
        rmsd = superimposer.rms
        print(f"Alignment RMSD: {rmsd:.3f} Å using {len(pdb_atoms)} CA atoms")
        
        return superimposer, rmsd
    
    def copy_residue(self, source_residue, target_chain):
        """
        Copy a residue from source to target chain.
        
        Args:
            source_residue: Source residue object
            target_chain: Target chain object
        """
        # Create new residue
        new_residue = Residue.Residue(
            source_residue.id,
            source_residue.resname,
            source_residue.segid
        )
        
        # Copy all atoms
        for atom in source_residue:
            new_atom = Atom.Atom(
                atom.name,
                atom.coord.copy(),
                atom.bfactor,
                atom.occupancy,
                atom.altloc,
                atom.fullname,
                atom.serial_number,
                atom.element
            )
            new_residue.add(new_atom)
        
        # Add residue to chain
        target_chain.add(new_residue)
    
    def complete_structure(self, chain_mapping=None):
        """
        Complete the PDB structure using AlphaFold predictions.
        
        Args:
            chain_mapping (dict): Mapping between PDB and AlphaFold chain IDs
                                 If None, assumes same chain IDs
        """
        if not self.pdb_structure or not self.af_structure:
            print("Error: Structures not loaded")
            return
        
        # Create a new structure for output
        completed_structure = Structure.Structure("completed")
        completed_model = Model.Model(0)
        completed_structure.add(completed_model)
        
        pdb_model = self.pdb_structure[0]
        af_model = self.af_structure[0]
        
        total_added = 0
        
        for pdb_chain in pdb_model:
            chain_id = pdb_chain.id
            
            # Find corresponding AlphaFold chain
            if chain_mapping and chain_id in chain_mapping:
                af_chain_id = chain_mapping[chain_id]
            else:
                af_chain_id = chain_id
            
            if af_chain_id not in af_model:
                print(f"Warning: Chain {af_chain_id} not found in AlphaFold structure")
                # Copy PDB chain as-is
                completed_model.add(pdb_chain.copy())
                continue
            
            af_chain = af_model[af_chain_id]
            
            # Align structures
            superimposer, rmsd = self.align_structures(pdb_chain, af_chain)
            
            if rmsd > 5.0:  # Threshold for acceptable alignment
                print(f"Warning: Poor alignment for chain {chain_id} (RMSD: {rmsd:.3f} Å)")
            
            # Find missing residues
            missing_residues = self.find_missing_residues(pdb_chain, af_chain)
            print(f"Chain {chain_id}: Found {len(missing_residues)} missing residues")
            
            # Create new chain for completed structure
            completed_chain = Chain.Chain(chain_id)
            
            # Copy existing PDB residues
            for residue in pdb_chain:
                if residue.id[0] == ' ':  # Standard residue
                    self.copy_residue(residue, completed_chain)
            
            # Add missing residues from AlphaFold
            added_count = 0
            for res_num in missing_residues:
                try:
                    af_residue = af_chain[res_num]
                    
                    # Check B-factor/confidence if available
                    avg_bfactor = np.mean([atom.bfactor for atom in af_residue])
                    
                    # Only add high-confidence residues (B-factor > 50 for AlphaFold)
                    if avg_bfactor > 50:
                        self.copy_residue(af_residue, completed_chain)
                        added_count += 1
                    else:
                        print(f"Skipping low-confidence residue {res_num} (confidence: {avg_bfactor:.1f})")
                        
                except KeyError:
                    print(f"Warning: Residue {res_num} not found in AlphaFold chain")
                    continue
            
            print(f"Added {added_count} high-confidence residues to chain {chain_id}")
            total_added += added_count
            
            completed_model.add(completed_chain)
        
        # Save completed structure
        io = PDBIO()
        io.set_structure(completed_structure)
        io.save(self.output_file)
        
        print(f"\nCompleted structure saved to: {self.output_file}")
        print(f"Total residues added: {total_added}")
    
    def analyze_completion(self):
        """Analyze the completion results."""
        if not os.path.exists(self.output_file):
            print("Error: Completed structure file not found")
            return
        
        # Load completed structure
        completed_structure = self.parser.get_structure("completed", self.output_file)
        
        print("\n=== Completion Analysis ===")
        
        for model in completed_structure:
            for chain in model:
                residue_count = len([r for r in chain if r.id[0] == ' '])
                print(f"Chain {chain.id}: {residue_count} residues")
        
        # Compare with original PDB
        print("\nOriginal PDB residue counts:")
        for model in self.pdb_structure:
            for chain in model:
                residue_count = len([r for r in chain if r.id[0] == ' '])
                print(f"Chain {chain.id}: {residue_count} residues")


def main():
    parser = argparse.ArgumentParser(
        description="Complete PDB structures using AlphaFold predictions"
    )
    parser.add_argument("pdb_file", help="Input PDB structure file")
    parser.add_argument("alphafold_file", help="AlphaFold structure file")
    parser.add_argument("-o", "--output", help="Output file name")
    parser.add_argument("--chain-mapping", help="Chain mapping file (format: pdb_chain:af_chain per line)")
    parser.add_argument("--confidence-threshold", type=float, default=50.0,
                       help="Minimum confidence score for adding residues (default: 50.0)")
    
    args = parser.parse_args()
    
    # Parse chain mapping if provided
    chain_mapping = None
    if args.chain_mapping:
        chain_mapping = {}
        with open(args.chain_mapping, 'r') as f:
            for line in f:
                if ':' in line:
                    pdb_chain, af_chain = line.strip().split(':')
                    chain_mapping[pdb_chain] = af_chain
    
    # Create completer and run
    completer = PDBCompleter(args.pdb_file, args.alphafold_file, args.output)
    completer.load_structures()
    completer.complete_structure(chain_mapping)
    completer.analyze_completion()

if __name__ == "__main__":
    main()