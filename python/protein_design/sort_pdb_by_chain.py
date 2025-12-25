
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
import os


def sort_pdb_by_chain(input_file, output_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_file)
    
    # Create a new structure to hold sorted chains
    new_structure = Structure('sorted_protein')
    new_model = Model(0)
    new_structure.add(new_model)
    
    # Get all chains and sort by chain ID
    chains = list(structure[0].get_chains())
    chains.sort(key=lambda x: x.id)
    
    # Add sorted chains to new structure
    for chain in chains:
        new_model.add(chain.copy())
    
    # Write sorted structure
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(output_file)


def sort_pdb_by_chain_from_dir(pdb_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    
    for pdb_file in pdb_files:
        input_file = os.path.join(pdb_dir, pdb_file)
        output_file = os.path.join(output_dir, pdb_file)
        sort_pdb_by_chain(input_file, output_file)


if __name__ == "__main__":

    import argparse
    
    parser = argparse.ArgumentParser(description='Sort PDB files by chain ID.')
    parser.add_argument('input', type=str, help='Input PDB file or directory containing PDB files.')
    parser.add_argument('output', type=str, help='Output file or directory to save sorted PDB files.')
    
    args = parser.parse_args()
    
    if os.path.isdir(args.input):
        sort_pdb_by_chain_from_dir(args.input, args.output)
    else:
        sort_pdb_by_chain(args.input, args.output)

    