#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 00:11:21 2025

@author: hao

conda env: unidock
required: rdkit, meeko, pandas
"""

import os
from subprocess import run
from rdkit import Chem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
import pandas as pd
from tqdm import tqdm

# === CONFIG ===
os.chdir('/home/hao/docking_tools/VS_HPK1')

receptor_pdbqt = "7M0M.pdbqt"

output_dir = "vina_result"
log_file = "docking_log.txt"
result_csv = "docking_results.csv"

library_sdf = "/home/hao/docking_tools/cmpd_library/VS_Campaign/diverse_20k.sdf"
ligand_dir = "/home/hao/docking_tools/cmpd_library/VS_Campaign/diverse_20k_pdbqt"


center = (-2.93, -6.681, -9.034)  # docking box center (x, y, z)
size = (15, 15, 15) # docking box size in Å

os.makedirs(ligand_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)




# === STEP 1: Convert SDF to individual PDBQT files ===
print("Converting SDF to PDBQT ligands...")
# supplier = Chem.SDMolSupplier(library_sdf)
# prep = MoleculePreparation()

# ligand_paths = []

# for i, mol in enumerate(supplier):
#     if mol is None:
#         continue
#     AllChem.EmbedMolecule(mol)
#     name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"ligand_{i+1}"
#     name = name.replace(" ", "_")
#     sdf_block = Chem.MolToMolBlock(mol)
#     prep.prepare(sdf_block)
#     out_path = os.path.join(ligand_dir, f"{name}.pdbqt")
#     with open(out_path, "w") as f:
#         f.write(prep.write_pdbqt_string())
#     ligand_paths.append((name, out_path))
    


supplier = Chem.SDMolSupplier(library_sdf, removeHs=False)
# iterate over molecules in SD file
for mol in tqdm(supplier, total=len(supplier)):
    pass
    mk_prep = MoleculePreparation()
    molsetup_list = mk_prep(Chem.AddHs(mol))
    molsetup = molsetup_list[0]
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)
    



# === STEP 3: Run Uni-Dock for each ligand ===
print("Running Uni-Dock...")
with open(log_file, "w") as logf:
    for name, ligand_path in ligand_paths:
        out_path = os.path.join(output_dir, f"{name}_out.pdbqt")
        log_path = os.path.join(output_dir, f"{name}_log.txt")

        cmd = [
            "unidock",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_path,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--output", out_path,
            "--log", log_path,
        ]

        run(cmd, stdout=logf, stderr=logf)

# === STEP 4: Parse log files and compile results ===
print("Parsing logs...")
results = []

for name, _ in ligand_paths:
    log_path = os.path.join(output_dir, f"{name}_log.txt")
    if not os.path.exists(log_path):
        continue

    with open(log_path) as f:
        for line in f:
            if "REMARK VINA RESULT" in line:
                parts = line.strip().split()
                if len(parts) >= 4:
                    score = float(parts[3])
                    results.append({"Ligand": name, "Affinity (kcal/mol)": score})
                break

df = pd.DataFrame(results)
df.sort_values(by="Affinity (kcal/mol)", inplace=True)
df.to_csv(result_csv, index=False)

print(f"Done! Results saved to {result_csv}")
