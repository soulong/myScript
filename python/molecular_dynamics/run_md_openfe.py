
# env: openfe

import sys, os, pathlib
import gufe
from openfe import NonTransformation, ChemicalSystem, ProteinComponent, SmallMoleculeComponent, SolventComponent
from openff.units import unit
from openfe.protocols.openmm_md.plain_md_methods import PlainMDProtocol

sys.path.append(os.path.expanduser('~/Documents/GitHub/myScript/python_functions/'))
from md_helper import prepare_protein, prepare_ligand_from_sdf, fix_trajectory


os.chdir('/media/hao/Data/Project_MYC/2025-07-18_AURKA_MYC_MD')

# prepare protein and ligand
if False:
    prepare_protein("5G1X_AURKA_NMYC.pdb", save_file="AURKA_NMYC.pdb", ignore_nonstandard_residues=True)
    prepare_protein("5G1X_AURKA.pdb", save_file="AURKA.pdb", ignore_nonstandard_residues=True)
    prepare_ligand_from_sdf('5G1X_OTS_diffdock.sdf', name="OTS", 
                            smiles="CC(=O)C1=CN=C2C=CC(=NC2=C1NC3CCC(CC3)CN(C)C)C4=CC(=C(C(=C4)Cl)O)Cl",
                            save_file='OTS.sdf')



# directory for md
md_dir = pathlib.Path("AURKA_NMYC")
md_dir.mkdir(parents=True, exist_ok=True)

protein_path = 'AURKA_NMYC.pdb'
# ligand_path = 'OTS.sdf'


protein = ProteinComponent.from_pdb_file(protein_path)
solvent = SolventComponent(ion_concentration=0.15 * unit.molar)

system = ChemicalSystem({'protein': protein, 'solvent': solvent})
# ligand = SmallMoleculeComponent.from_sdf_file(ligand_path)
# system = ChemicalSystem({'ligand': ligand, 'protein': protein, 'solvent': solvent})


# set up the protocol
settings = PlainMDProtocol.default_settings()
settings.protocol_repeats = 1
settings.engine_settings.compute_platform = 'cuda' # Use the default best platform
settings.solvation_settings.solvent_padding = 1.0 * unit.nanometer
settings.integrator_settings.timestep = 4 * unit.femtosecond
settings.simulation_settings.minimization_steps = 50000
settings.simulation_settings.equilibration_length_nvt = 1 * unit.nanosecond
settings.simulation_settings.equilibration_length = 1 * unit.nanosecond 
settings.simulation_settings.production_length = 300 * unit.nanosecond
settings.output_settings.checkpoint_interval = 1000 * unit.picosecond
settings.output_settings.trajectory_write_interval = 40 * unit.picosecond
settings.output_settings.output_indices = 'protein or resname UNK'
print(settings)
protocol = PlainMDProtocol(settings=settings)


## option 1
nontransformation = NonTransformation(system=system, protocol=protocol, name=f"{system.name}")
nontransformation.to_json(md_dir / f"protocol.json")
bash_command = f'openfe quickrun {str(md_dir)}/protocol.json -o {str(md_dir)}/results.json -d {str(md_dir)}'
print(bash_command) # run in terminal, conda activate openfe, cd to the directory

# ## option 2
# dag = protocol.create(stateA=system, stateB=system, mapping=None)
# dagres = gufe.protocols.execute_DAG(
#     dag,
#     shared_basedir=md_dir, scratch_basedir=md_dir,
#     keep_shared=True, # set this to True to save the outputs
#     n_retries=1
#     )




# Post-process trajectory
if False:
    import os
    os.chdir('')
    fix_trajectory(
        top='complex.pdb',
        traj='simulation.xtc',                  
        atomgroup="not (resname HOH or resname SOL or resname NA or resname CL)", 
        out_traj='simulation_clean.xtc',
        out_pdb='frame.pdb'
        )




