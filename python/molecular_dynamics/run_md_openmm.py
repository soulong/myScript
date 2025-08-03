
# env: openmm

import os, sys
from openmm import Platform, MonteCarloBarostat
from openmm import app, unit, LangevinMiddleIntegrator
from openmm.app import PDBFile, Simulation, Modeller, StateDataReporter, XTCReporter, CheckpointReporter
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator

sys.path.append(os.path.expanduser('~/Documents/GitHub/myScript/python_functions/'))
from md_helper import prepare_protein, prepare_ligand_from_sdf, clean_trajectory

os.chdir('/media/hao/Data/Project_MYC/2025-07-18_AURKA_MYC_MD')



# prepare_protein("5G1X_clean.pdb", save_file="protein_prepared.pdb")
# prepare_ligand_from_sdf('OTSSP167_diffdock.sdf', "OTS", 
#                         "CC(=O)C1=CN=C2C=CC(=NC2=C1NC3CCC(CC3)CN(C)C)C4=CC(=C(C(=C4)Cl)O)Cl",
#                         save_file='ligand_prepared.sdf')


# directory for md
output_dir = "AURKA_NMYC_OTS"

protein_path = 'protein_prepared.pdb'
ligand_path = 'ligand_prepared.sdf'
# ligand_path = None


os.makedirs(output_dir, exist_ok=True)
os.system(f"cp {protein_path} {output_dir}/{protein_path}")
if ligand_path: os.system(f"cp {ligand_path} {output_dir}/{ligand_path}")
protein_path = os.path.join(output_dir, protein_path)
ligand_path = os.path.join(output_dir, ligand_path) if ligand_path else None

protein = PDBFile(protein_path)
ligand = Molecule.from_file(ligand_path) if ligand_path else None


# set up SystemGenerator
# SystemGenerator().SMALL_MOLECULE_FORCEFIELDS
system_generator = SystemGenerator(
    forcefields = ['amber/ff14SB.xml', 'amber/phosaa14SB.xml','amber/tip3p_standard.xml'],
    small_molecule_forcefield = 'openff-2.2.0', # gaff-2.11
    molecules = [ligand, ],
    forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 
                        'ewaldErrorTolerance': 0.0005, 'hydrogenMass': 3*unit.amu},
    periodic_forcefield_kwargs = {'nonbondedMethod': app.PME,
                                  'nonbondedCutoff': 1 * unit.nanometer},
    cache=f'{output_dir}/db.json')

# Use Modeller to combine the protein and ligand into a complex
modeller = Modeller(protein.topology, protein.positions)
if ligand:
    ligand_top = ligand.to_topology()
    modeller.add(ligand_top.to_openmm(), ligand_top.get_positions().to_openmm())
    print(modeller.getTopology())
    
modeller.addSolvent(system_generator.forcefield, model='tip3p', 
                    boxShape='cube', padding=10 * unit.angstroms,
                    positiveIon='Na+', negativeIon='Cl-', ionicStrength=0.15 * unit.molar)
print(modeller.getTopology())

# save the complex
output_complex=f'{output_dir}/complex.pdb'
with open(output_complex, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# set up system
system = system_generator.create_system(modeller.topology, molecules=[ligand, ] if ligand else None)
system.addForce(MonteCarloBarostat(1 * unit.atmospheres, 300 * unit.kelvin, 25))

# set up integrator
step_size = 0.004 # ps
temperature = 300  # kelvin
integrator = LangevinMiddleIntegrator(temperature * unit.kelvin, 1 / unit.picosecond, step_size * unit.picoseconds)
integrator.setConstraintTolerance(0.000001)

# set up simulation
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(modeller.positions)


# minimize
print('Minimising ...')
simulation.minimizeEnergy()
with open(f'{output_dir}/minimized.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)

# equilibrate
print('Equilibrating ...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(250000)  # 250000 steps at 0.004 ps per step = 1000 ps = 1 ns
equilibrated=f'{output_dir}/equilibrated.pdb'
with open(equilibrated, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)

# simulation
duration = 200 # ns
reporting_interval = int((duration * 1000 / step_size) / 5000)  # total reports = 5000
num_steps = int(duration * 1000 / step_size)  # total steps = duration in ps / step size in ps
output_traj = f'{output_dir}/trajectory.xtc'
simulation.reporters.append(
    XTCReporter(output_traj, reporting_interval, enforcePeriodicBox=True))
simulation.reporters.append(
    CheckpointReporter(f'{output_dir}/checkpoint20.chk', int(reporting_interval / 20)))
simulation.reporters.append(
    StateDataReporter(f'{output_dir}/log.txt', reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(
    StateDataReporter(sys.stdout, reporting_interval * 10, step=True, potentialEnergy=True, temperature=True))
print('Starting simulation with', num_steps, 'steps ...')
simulation.currentStep = 0
simulation.step(num_steps)
simulation.saveState(f'{output_dir}/final_state.xml')


# clean up
clean_trajectory(output_traj, equilibrated,
                 output_traj=f'{output_dir}/trajectory_clean.xtc', 
                 output_complex=f'{output_dir}/complex_clean.pdb')