
# import os
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import copy
from IPython.display import display
import pdbfixer
from openmm.app import PDBFile

# import mdtraj as md
import MDAnalysis as mda
import MDAnalysis.transformations as trans


def prepare_protein(pdb_path, ignore_nonstandard_residues=False, ph=7.0, save_file=None):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """

    fixer = pdbfixer.PDBFixer(pdb_path)
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer

    if ignore_nonstandard_residues:
        fixer.findNonstandardResidues()  # find non-standard residue
        fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one

    fixer.findMissingResidues() # identify missing residues, needed for identification of missing atoms
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues

    fixer.addMissingHydrogens(ph)  # add missing hydrogens

    if save_file:
        with open(save_file, 'w') as outfile:
            PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile, keepIds=True)

    return fixer

def prepare_ligand_from_pdb(pdb_path, resname, name=None, smiles=None, save_file=None, depict=False):
    """
    Prepare a ligand from a PDB file via adding hydrogens and assigning bond orders.

    Parameters
    ----------
    pdb_path: str
       PDB file containing the ligand of interest.
    resname: str
        Three character residue name of the ligand.
    smiles : str
        SMILES string of the ligand informing about correct protonation and bond orders.
    depict: bool, optional
        show a 2D representation of the ligand

    Returns
    -------
    prepared_ligand: rdkit.Chem.Mol
        Prepared ligand.
    """

    # split molecule
    rdkit_mol = Chem.MolFromPDBFile(pdb_path, sanitize=True, removeHs=True)
    rdkit_mol_split = Chem.rdmolops.SplitMolByPDBResidues(rdkit_mol)

    # extract the ligand and remove any already present hydrogens
    ligand = rdkit_mol_split[resname]
    ligand = Chem.RemoveHs(ligand)

    # assign bond orders from template
    if smiles:
        # if a SMILES string is provided, use it as a template for bond orders
        reference_mol = Chem.MolFromSmiles(smiles)
        prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference_mol, ligand)
        # prepared_ligand.AddConformer(ligand.GetConformer(0))
    else:
        prepared_ligand = ligand

    if name is None:
        prepared_ligand.SetProp("_Name", name)

    prepared_ligand = Chem.AddHs(prepared_ligand, addCoords=True)

    if save_file:
        with Chem.SDWriter(save_file) as w:
            w.write(prepared_ligand)
    
    if depict:
        ligand_2d = copy.deepcopy(ligand)
        prepared_ligand_2d = copy.deepcopy(prepared_ligand)
        AllChem.Compute2DCoords(ligand_2d)
        AllChem.Compute2DCoords(prepared_ligand_2d)
        display(
            Draw.MolsToGridImage(
                [ligand_2d, prepared_ligand_2d], molsPerRow=2, legends=["original", "prepared"]
            )
        )

    return prepared_ligand

def prepare_ligand_from_sdf(sdf_path, name=None, smiles=None, save_file=None, depict=False):
    """
    Prepare a ligand from an SDF file via adding hydrogens and assigning bond orders.

    Parameters
    ----------
    sdf_path: str
        SDF file containing the ligand of interest.
    name: str, optional
        Name of the ligand to be set in the SDF file.
    smiles : str, optional
        SMILES string of the ligand informing about correct protonation and bond orders.
    save_file: str, optional
        Path to save the prepared ligand in SDF format.
    depict: bool, optional
        show a 2D representation of the ligand
    
    Returns
    -------
    prepared_ligand: rdkit.Chem.Mol
        Prepared ligand.
    """

    ligand = Chem.MolFromMolFile(sdf_path, sanitize=True, removeHs=True)

    if smiles:
        reference_mol = Chem.MolFromSmiles(smiles,sanitize=True)
        prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference_mol, ligand)
        # prepared_ligand.AddConformer(ligand.GetConformer(0))
    else:
        prepared_ligand = ligand

    if name is None:
        prepared_ligand.SetProp("_Name", name)

    prepared_ligand = Chem.AddHs(prepared_ligand, addCoords=True)

    if save_file:
        with Chem.SDWriter(save_file) as w:
            w.write(prepared_ligand)

    if depict:
        ligand_2d = copy.deepcopy(ligand)
        prepared_ligand_2d = copy.deepcopy(prepared_ligand)
        AllChem.Compute2DCoords(ligand_2d)
        AllChem.Compute2DCoords(prepared_ligand_2d)
        display(
            Draw.MolsToGridImage(
                [ligand_2d, prepared_ligand_2d], molsPerRow=2, legends=["original", "prepared"]
            )
        )

    return prepared_ligand




# Load the trajectory and topology
# os.chdir('/media/hao/Data/Project_MYC/2025-07-18_AURKA_MYC_MD/AURKA_NMYC_OTS/shared_PlainMDProtocolUnit-b8183f09785a4938bf787a9a6ee70ec8_attempt_0')
def fix_trajectory(
    top: str = "equil_npt.pdb",
    traj: str = "simulation.xtc",
    atomgroup: str = "not (resname HOH or resname SOL or resname NA or resname CL)",
    out_traj: str = "simulation_clean.xtc",
    out_pdb: str = "frame.pdb",
    thread: int = 10
    ) -> mda.Universe:
    """
    Fixes periodic boundary conditions and rotational alignment for a molecular dynamics trajectory.
    Removes solvent and ions, unwraps molecules, centers and aligns the system, and writes cleaned outputs.

    Parameters
    ----------
    top : str
    Path to the topology file (e.g., PDB).
    traj : str
    Path to the trajectory file (e.g., XTC).
    atomgroup : str
    Atom selection string for MDAnalysis to specify which atoms to keep (e.g., exclude water and ions).
    out_traj : str
    Path to save the cleaned trajectory file.
    out_pdb : str
    Path to save a PDB file of the selected atoms (first frame).
    thread : int
    Number of threads to use for transformations.

    Returns
    -------
    u : MDAnalysis.Universe
    The MDAnalysis Universe object with applied transformations.
    """

    print("read trajectory")
    u = mda.Universe(top, traj, in_memory=True)

    print("chains:", *u.segments.segids)

    # Select the atoms to unwrap (e.g., all atoms)
    ag = u.select_atoms(atomgroup)

    # Add the nojump transformation
    print("apply transformation")
    u.trajectory.add_transformations(
    trans.NoJump(max_threads=thread), 
    trans.center_in_box(ag, max_threads=thread), 
    trans.fit_rot_trans(ag, ag, weights="mass", max_threads=thread)
    )

    if out_traj:
        print("write updated trajectory")
        with mda.Writer(out_traj, multiframe=True) as W:
            # Iterate through the trajectory
            for ts in u.trajectory:
                # Write the current frame to the file
                W.write(ag)

    if out_pdb:
        ag.write(out_pdb)

    return u







import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms, pca
import prolif as plf
from rdkit import DataStructs
from matplotlib import pyplot as plt
import seaborn as sns
import nglview as nv
import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any


class MDAnalysisHelper:
    """A helper class for performing molecular dynamics analysis using MDAnalysis."""
    
    def __init__(self, topology: str = "frame.pdb", 
                 trajectory: str = "simulation_clean.xtc", 
                 in_memory: bool = True):
        """
        Initialize the MDAnalysisHelper with a topology and trajectory file.
        
        Args:
            topology (str): Path to the topology file (default: "frame.pdb").
            trajectory (str): Path to the trajectory file (default: "simulation_clean.xtc").
            in_memory (bool): Whether to load the trajectory into memory for faster processing (default: True).
        
        Raises:
            FileNotFoundError: If the topology or trajectory file does not exist.
            ValueError: If the files are not compatible with MDAnalysis.
        """
        try:
            self.universe = mda.Universe(topology, trajectory, in_memory=in_memory)
            self.chains = self.universe.segments.segids.tolist()
        except FileNotFoundError as e:
            raise FileNotFoundError(f"Could not find file: {e}")
        except Exception as e:
            raise ValueError(f"Error loading files with MDAnalysis: {e}")

    def run_rmsd(self,
                 select: str = 'backbone',
                 ag_str_list: List[str] = ['protein'],
                 name_list: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Calculate the Root Mean Square Deviation (RMSD) over a trajectory.
        
        Args:
            select (str): Selection string for alignment (default: 'backbone').
            ag_str_list (List[str]): List of selection strings for RMSD calculation (default: ['protein']).
            name_list (Optional[List[str]]): Custom names for the selections in output (default: None).
        
        Returns:
            pd.DataFrame: DataFrame containing time and RMSD values for each selection.
        
        Raises:
            ValueError: If the selection strings are invalid.
        """
        try:
            self.universe.trajectory[0]  # Set to first frame
            analysis = rms.RMSD(self.universe, select=select, groupselections=ag_str_list).run()
            name_list = ag_str_list if name_list is None else name_list
            result = pd.DataFrame(analysis.results.rmsd, columns=['step', 'time', select, *name_list])
            return result.drop('step', axis=1)
        except Exception as e:
            raise ValueError(f"Error in RMSD calculation: {e}")

    def run_rmsf(self, ag_str: str = 'protein and name CA') -> pd.DataFrame:
        """
        Calculate the Root Mean Square Fluctuation (RMSF) for selected atoms.
        
        Args:
            ag_str (str): Selection string for RMSF calculation (default: 'protein and name CA').
        
        Returns:
            pd.DataFrame: DataFrame with residue information and RMSF values.
        
        Raises:
            ValueError: If the selection string is invalid or RMSF calculation fails.
        """
        try:
            ag = self.universe.select_atoms(ag_str)
            if len(ag) == 0:
                raise ValueError(f"No atoms selected with: {ag_str}")
            resid = ag.residues.resids.tolist()
            resname = ag.residues.resnames.tolist()
            segid = ag.residues.segids.tolist()
            rmsf = rms.RMSF(ag).run().results.rmsf.tolist()
            result = pd.DataFrame({'resid': resid, 'resname': resname, 'segid': segid, 'rmsf': rmsf})
            
            # Set B-factor for visualization
            self.universe.add_TopologyAttr('tempfactors')
            protein = self.universe.select_atoms('protein')
            for residue, r_value in zip(protein.residues, rmsf):
                residue.atoms.tempfactors = r_value
                
            return result
        except Exception as e:
            raise ValueError(f"Error in RMSF calculation: {e}")

    def run_prolif(self,
                   ligand: str = 'segid A',
                   protein: str = 'protein and byres around 20.0 group ligand') -> Dict[str, Any]:
        """
        Run Protein-Ligand Interaction Fingerprint (ProLIF) analysis.
        
        Args:
            ligand (str): Selection string for ligand (default: 'segid A').
            protein (str): Selection string for protein (default: 'protein and byres around 20.0 group ligand').
        
        Returns:
            Dict[str, Any]: Dictionary containing fingerprint object, DataFrame, and similarity matrix.
        
        Raises:
            ValueError: If selections are invalid or ProLIF analysis fails.
        """
        try:
            ligand_selection = self.universe.select_atoms(ligand)
            if len(ligand_selection) == 0:
                raise ValueError(f"No ligand atoms selected with: {ligand}")
            protein_selection = self.universe.select_atoms(protein, ligand=ligand_selection)
            if len(protein_selection) == 0:
                raise ValueError(f"No protein atoms selected with: {protein}")
                
            fp = plf.Fingerprint().run(self.universe.trajectory, ligand_selection, protein_selection)
            bitvectors = fp.to_bitvectors()
            df = fp.to_dataframe()
            similarity_matrix = pd.DataFrame(
                [DataStructs.BulkTanimotoSimilarity(bv, bitvectors) for bv in bitvectors],
                index=df.index, columns=df.index
            )
            return {'fp': fp, 'df': df, 'similarity_matrix': similarity_matrix}
        except Exception as e:
            raise ValueError(f"Error in ProLIF analysis: {e}")

    def run_pca(self, select_str: str = 'backbone', n_components: int = 2) -> Dict[str, Any]:
        """
        Perform Principal Component Analysis (PCA) on the trajectory.
        
        Args:
            select_str (str): Selection string for PCA (default: 'backbone').
            n_components (int): Number of principal components to compute (default: 2).
        
        Returns:
            Dict[str, Any]: Dictionary containing PCA object and transformed DataFrame.
        
        Raises:
            ValueError: If PCA calculation fails or selection is invalid.
        """
        try:
            pc = pca.PCA(self.universe, select=select_str, align=True, mean=None, n_components=n_components).run()
            for i in range(n_components):
                print(f"Cumulated variance: {pc.cumulated_variance[i]:.3f}")
            transformed = pc.transform(self.universe.select_atoms(select_str), n_components=n_components)
            df = pd.DataFrame(transformed, columns=[f'PC{i+1}' for i in range(n_components)])
            df['time'] = df.index * self.universe.trajectory.dt
            return {'pc': pc, 'df': df}
        except Exception as e:
            raise ValueError(f"Error in PCA calculation: {e}")





if __name__ == "main":
    
    root_dir='/media/hao/Data/Project_MYC/2025-07-18_AURKA_MYC_MD/'
    os.chdir(f'{root_dir}/AURKA_NMYC')
    os.getcwd()

    u = fix_trajectory(
        top = "equil_npt.pdb",
        traj = "simulation.xtc",
        atomgroup = "not (resname HOH or resname SOL or resname NA or resname CL)"
        )

    ana = MDAnalysisHelper('frame.pdb', 'simulation_clean.xtc', in_memory=True)

    # rmsd
    ana_rmsd = ana.run_rmsd(ag_str_list=[f'segid {x}' for x in ana.universe.segments.segids.tolist()])
    ana_rmsd.to_csv('rmsd.csv', index=True)

    ana_rmsd = ana_rmsd.set_index('time')
    ana_rmsd.plot(xlabel='Time (ps)', ylabel="RMSD (A")
    # plt.show()
    fig = plt.gcf()
    fig.set_size_inches(8, 5)
    plt.savefig("rmsd.pdf")

    # rmsf
    ana_rmsf = ana.run_rmsf()
    ana_rmsf.to_csv('rmsf.csv', index=True)

    ana_rmsf['uid'] = ana_rmsf['segid'].astype(str) + '.' + ana_rmsf['resid'].astype(str) + '.' + ana_rmsf['resname'].astype(str) 
    # ana_rmsf = ana_rmsf.set_index('uid')
    ax = ana_rmsf.plot(x='uid', y='rmsf', xlabel='', ylabel="RMSD (A")
    plt.xticks(range(len(ana_rmsf['uid'])), ana_rmsf['uid'], rotation=90, fontsize=4, ha='right')
    fig = plt.gcf()
    fig.set_size_inches(12, 5)
    plt.savefig("rmsf.pdf")

    # prolif
    # A-B
    prolif_AB = ana.run_prolif(
        ligand='segid B', protein='protein and byres around 20.0 group ligand')    
    prolif_AB['df'].to_csv('prolif_AB.csv')

    prolif_AB['fp'].plot_barcode()
    fig = plt.gcf()
    fig.set_size_inches(10, 15)
    plt.savefig("prolif_AB.pdf")

    # A-X
    prolif_AX = ana.run_prolif(
        ligand='segid X', protein='protein and byres around 20.0 group ligand')    
    prolif_AX['df'].to_csv('prolif_AX.csv')

    prolif_AX['fp'].plot_barcode()
    fig = plt.gcf()
    fig.set_size_inches(10, 10)
    plt.savefig("prolif_AX.pdf")

    # pca
    pc = ana.run_pca()
