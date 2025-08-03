import os
from pathlib import Path
from typing import Dict
from Bio import SeqIO
import pandas as pd

import pyrosetta
from pyrosetta import pose_from_pdb, standard_packer_task
from pyrosetta.rosetta.core.pose import remove_nonprotein_residues
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.pack.task.operation import (
    RestrictToRepacking, PreventRepackingRLT, OperateOnResidueSubset, InitializeFromCommandline
)
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover, MinMover
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.core.pose import get_chain_from_chain_id
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue


def mutate_pdb_with_sequence(pdb_file, sequence, output_file, minimize=False, chain_id="A"):
    # pyrosetta.init("-mute all")
    pose = pose_from_pdb(pdb_file)

    # Clean non-protein atoms (e.g., water, ligands)
    remove_nonprotein_residues(pose)

    # Select chain residues
    selector = ChainSelector(chain_id)
    selected_residues = [i + 1 for i, sel in enumerate(selector.apply(pose)) if sel]

    if len(sequence) != len(selected_residues):
        raise ValueError(f"Sequence length ({len(sequence)}) doesn't match number of residues in chain {chain_id} ({len(selected_residues)}).")

    print(f"Mutating {len(selected_residues)} residues in chain {chain_id}...")

    # Perform mutations
    # for i, aa in zip(selected_residues, sequence):
    #     target_aa = aa_from_oneletter_code(aa.upper())
    #     pose.replace_residue(i, pyrosetta.rosetta.core.conformation.ResidueFactory.create_residue(
    #         pose.residue(i).residue_type_set().name_map(target_aa)), True)
    for i, aa in zip(selected_residues, sequence):
        print(f"Residue {i}: {pose.residue(i).name1()}, Chain {pose.pdb_info().chain(i)}")
        target_aa = aa.upper()
        if pose.residue(i).name1() == target_aa:
            continue
        mutator = MutateResidue(i, target_aa)
        mutator.apply(pose)

    # Packing task
    print("packing sidechain")
    tf = standard_packer_task(pose)
    tf.restrict_to_repacking()
    tf.or_include_current(True)

    # Prevent repacking other chains
    for cid in set(pose.pdb_info().chain(i) for i in range(1, pose.size() + 1)):
        if cid != chain_id:
            tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(), ChainSelector(cid), True))

    # Sidechain packing
    pack_mover = PackRotamersMover()
    pack_mover.task_factory(tf)
    pack_mover.score_function(pyrosetta.get_fa_scorefxn())
    pack_mover.apply(pose)

    # Optional: Minimize sidechains
    if minimize:
        print(f"Minimizing sidechains for chain {chain_id}...")
        mm = MoveMap()
        mm.set_bb(False)
        mm.set_chi(False)
        for res in selected_residues:
            mm.set_chi(res, True)

        scorefxn = pyrosetta.get_fa_scorefxn()
        min_mover = MinMover()
        min_mover.movemap(mm)
        min_mover.score_function(scorefxn)
        min_mover.min_type("dfpmin")  # or "lbfgs_armijo_nonmonotone"
        min_mover.tolerance(0.0001)
        min_mover.apply(pose)

    # Save mutated and repacked structure
    pose.dump_pdb(output_file)
    print(f"Output saved to: {output_file}")


def load_sequences(input_file: str) -> Dict[str, str]:
    sequences = {}
    if input_file.endswith(".csv"):
        df = pd.read_csv(input_file)
        if 'id' not in df.columns or 'sequence' not in df.columns:
            raise ValueError("CSV must contain 'id' and 'sequence' columns.")
        for _, row in df.iterrows():
            sequences[row['id']] = row['sequence']
    elif input_file.endswith(".fasta") or input_file.endswith(".fa"):
        for record in SeqIO.parse(input_file, "fasta"):
            sequences[record.id] = str(record.seq)
    else:
        raise ValueError("Input sequence file must be CSV or FASTA format.")
    return sequences


def batch_mutate_all(backbone_dir: str, sequence_file: str, output_dir: str, chain_id: str = "A", minimize=False):
    sequences = load_sequences(sequence_file)
    backbone_dir = Path(backbone_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for seq_id, seq in sequences.items():
        if "#" not in seq_id:
            print(f"Skipping malformed sequence ID (no '#'): {seq_id}")
            continue
        backbone_name, tag = seq_id.split("#", 1)
        pdb_path = backbone_dir / f"{backbone_name}.pdb"
        if not pdb_path.exists():
            print(f"[Warning] PDB file not found for {backbone_name}. Skipping...")
            continue

        output_path = output_dir / f"{backbone_name}_{tag}.pdb"
        print(f"Processing: {seq_id} -> {output_path.name} [Chain {chain_id}]")
        try:
            mutate_pdb_with_sequence(str(pdb_path), seq, str(output_path), minimize=minimize, chain_id=chain_id)
        except Exception as e:
            print(f"[Error] Failed to process {seq_id}: {e}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Batch mutate and repack multiple backbone PDBs with designed sequences.")
    parser.add_argument("backbone_dir", help="Directory containing backbone PDB files (e.g., backbone_0.pdb)")
    parser.add_argument("sequence_file", help="CSV or FASTA file containing sequences (ID format: backboneName#tag)")
    parser.add_argument("output_dir", help="Directory to write mutated output PDBs")
    parser.add_argument("--chain", default="A", help="Chain ID in backbone PDB to mutate (default: A)")
    parser.add_argument("--minimize", action="store_true", help="Minimize sidechains after repacking")

    args = parser.parse_args()
    pyrosetta.init("-mute all -out:level 0")
    batch_mutate_all(args.backbone_dir, args.sequence_file, args.output_dir, chain_id=args.chain, minimize=args.minimize)
