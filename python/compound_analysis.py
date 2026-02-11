## env: chem
## date: 2026-01-05

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit.SimDivFilters import MaxMinPicker
import umap
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg
from io import BytesIO 
import os
import joblib
from datetime import datetime
from typing import Union, Optional, List

# Suppress the UMAP Jaccard gradient warning to keep console clean
import warnings
warnings.filterwarnings("ignore", message="gradient function is not yet implemented")

class ChemicalLibraryAnalyzer:
    """
    Core class for the full library: processes once, computes UMAP and diversity clusters.
    Now includes high-speed Tanimoto similarity searching.
    """
    
    def __init__(
        self,
        id_col: str = 'compound_id',
        smiles_col: str = 'smiles'
    ):
        self.id_col = id_col
        self.smiles_col = smiles_col
        
        self.df = None              # Library dataframe
        self.fps_array = None       # For UMAP
        self.diverse_indices = None # Diversity representatives

    def __repr__(self):
        status = "Empty" if self.df is None else f"{len(self.df)} compounds"
        return f"<ChemicalLibraryAnalyzer: {status}>"
    
    def _load_from_file(self, file_path: str, ) -> pd.DataFrame:
        if file_path.endswith('.csv'):
            return pd.read_csv(file_path)
        elif file_path.endswith(('.xlsx', '.xls')):
            return pd.read_excel(file_path)
        else:
            raise ValueError("File must be .csv or .xlsx")
    
    def load_library(
        self,
        source: Union[str, pd.DataFrame],
        id_col: Optional[str] = None,
        smiles_col: Optional[str] = None
    ):
        if isinstance(source, pd.DataFrame):
            df = source.copy()
        else:
            df = self._load_from_file(source)
        
        if id_col:
            self.id_col = id_col
        if smiles_col:
            self.smiles_col = smiles_col
        
        df.columns = df.columns.str.strip()
        df = df.rename(columns={self.id_col: 'compound_id', self.smiles_col: 'smiles'})

        required = ['compound_id', 'smiles']
        missing = [c for c in required if c not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        
        self.df = df[required].astype({'compound_id': str})
        print(f"[*] Loaded library with {len(self.df)} compounds.")
    
    def generate_fingerprints_and_properties(
        self,
        radius: int = 2,
        n_bits: int = 2048
    ):
        print(f"[*] Parsing SMILES and calculating properties (Radius={radius})...")
        def process(smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            mw = Descriptors.MolWt(mol)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            return (mol, mw, fp)
        
        results = self.df['smiles'].apply(process)
        valid = results.notna()
        self.df = self.df[valid].reset_index(drop=True)
        
        processed = results[valid].tolist()
        self.df['mol'] = [x[0] for x in processed]
        self.df['mw'] = [x[1] for x in processed]
        self.df['fp'] = [x[2] for x in processed]
        
        print(f"[*] {len(self.df)} valid molecules remain.")
        
        # Fingerprint matrix for UMAP (Boolean for Jaccard efficiency)
        self.fps_array = np.zeros((len(self.df), n_bits), dtype=bool)
        for i, fp in enumerate(self.df['fp']):
            self.fps_array[i] = fp.ToList()
    def select_diverse_representatives(
        self,
        num_reps: int = 300,
        seed: int = 42
    ):
        print(f"[*] Selecting {num_reps} diverse representatives...")
        picker = MaxMinPicker()
        fps_list = self.df['fp'].tolist()
        self.diverse_indices = list(picker.LazyBitVectorPick(fps_list, len(fps_list), num_reps, seed=seed))
        
        # Assign every compound to nearest representative → cluster
        diverse_fps = [fps_list[i] for i in self.diverse_indices]
        clusters = []
        for fp in fps_list:
            sims = BulkTanimotoSimilarity(fp, diverse_fps)
            clusters.append(np.argmax(sims))
        
        self.df['cluster'] = clusters
        print(f"[*] Assigned all compounds to {num_reps} diversity clusters.")
    
    def run_umap(
        self,
        n_neighbors: int = 20,
        min_dist: float = 0.05,
        random_state: int = 42
    ):
        """
        Runs UMAP. Default parameters are optimized for cluster separation:
        Lower n_neighbors (20) focuses on local structure.
        Lower min_dist (0.05) allows tighter clumping.
        """
        print(f"[*] Running UMAP (n_neighbors={n_neighbors}, min_dist={min_dist})...")
        reducer = umap.UMAP(
            metric='jaccard',
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=2,
            random_state=random_state,
            low_memory=True,
            verbose=False
        )
        embedding = reducer.fit_transform(self.fps_array)
        self.df['umap_x'] = embedding[:, 0]
        self.df['umap_y'] = embedding[:, 1]
        print("[*] UMAP completed.")

    def search_similar(self, query_smiles: str, top_n: int = 10) -> pd.DataFrame:
        """
        Rank library compounds by Tanimoto similarity to a query SMILES.
        """
        query_mol = Chem.MolFromSmiles(query_smiles)
        if query_mol is None:
            raise ValueError(f"Invalid SMILES string: {query_smiles}")
        
        # Must use same params as library fingerprints
        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=2048)
        sims = BulkTanimotoSimilarity(query_fp, self.df['fp'].tolist())
        
        sim_series = pd.Series(sims, index=self.df.index)
        top_indices = sim_series.nlargest(top_n).index
        
        results = self.df.loc[top_indices].copy()
        results['tanimoto_similarity'] = sim_series.loc[top_indices]
        print(f"[*] Found {top_n} analogs. Max Similarity: {results['tanimoto_similarity'].max():.3f}")
        return results


class HitVisualizer:
    """
    Class to overlay hits on the pre-processed library.
    """
    
    def __init__(self, library_analyzer: ChemicalLibraryAnalyzer):
        self.lib = library_analyzer
        if 'umap_x' not in self.lib.df.columns:
            raise ValueError("Library must have UMAP coordinates (run_umap first).")
    
    def add_hits(
        self,
        source: Optional[Union[str, pd.DataFrame]] = None,
        id_col: Optional[str] = None,
        name: str = 'hits'
    ) -> pd.DataFrame:
        work_df = self.lib.df.copy()
        work_df['is_hit'] = False
        
        if source is None:
            return work_df
        
        if isinstance(source, pd.DataFrame):
            hits_df = source.copy()
        else:
            if source.endswith('.csv'):
                hits_df = pd.read_csv(source)
            elif source.endswith(('.xlsx', '.xls')):
                hits_df = pd.read_excel(source)
            else:
                raise ValueError("Hits file must be .csv or .xlsx")
        
        if id_col:
            hits_df = hits_df.rename(columns={id_col: 'compound_id'})
        
        hits_df['compound_id'] = hits_df['compound_id'].astype(str)
        merge_cols = [c for c in hits_df.columns if c != 'compound_id']
        
        if merge_cols:
            work_df = work_df.merge(hits_df[['compound_id'] + merge_cols], 
                                    on='compound_id', how='left')
        
        work_df['is_hit'] = work_df['compound_id'].isin(hits_df['compound_id'])
        print(f"[*] [{name}] Matched {work_df['is_hit'].sum()} hits.")
        return work_df
    
    def map_umap(
        self,
        work_df: pd.DataFrame,
        output_file: str = 'chemical_space.pdf',
        background_subsample: int = 25000,
        figsize: tuple = (12, 10),
        plot_title: Optional[str] = None,
        point_size: Optional[str] = None,           # NEW: column name to scale hit point size
        point_size_scale: float = 1.0,              # Multiplier for size range
        point_size_min: int = 20,                   # Minimum hit point size
        point_size_max: int = 150,                  # Maximum hit point size
        background_size: int = 8                    # Fixed size for background points
    ):
        """
        Enhanced UMAP plot with optional variable point sizing for hits.
        
        Parameters:
        - point_size: str or None
            Column name in work_df to use for scaling hit point sizes.
            Must be numeric. If None, all hits use fixed size (40).
        - point_size_scale: float
            Scaling factor to adjust overall size range.
        - point_size_min/max: int
            Bounds for scaled sizes.
        """
        with PdfPages(output_file) as pdf:
            fig, ax = plt.subplots(figsize=figsize)
            
            is_hit_present = work_df['is_hit'].any()
            
            # Background points (always small and fixed size)
            bg = work_df[~work_df['is_hit']].copy()
            if len(bg) > background_subsample:
                bg = bg.sample(n=background_subsample, random_state=42)
            
            if is_hit_present:
                sns.scatterplot(
                    data=bg, x='umap_x', y='umap_y',
                    color='lightgray', alpha=0.5, s=background_size,
                    ax=ax, rasterized=True, legend=False
                )
            else:
                # No hits: color by cluster
                sns.scatterplot(
                    data=bg, x='umap_x', y='umap_y', hue='cluster',
                    palette='tab20', alpha=0.6, s=background_size,
                    ax=ax, rasterized=True, legend=False
                )
            
            # Hits layer
            if is_hit_present:
                hits_df = work_df[work_df['is_hit']].copy()
                
                if point_size and point_size in hits_df.columns:
                    # Ensure numeric and handle NaN
                    size_col = hits_df[point_size].replace([np.inf, -np.inf], np.nan)
                    if size_col.dtype.kind not in 'bifc':  # not numeric
                        raise ValueError(f"point_size column '{point_size}' must be numeric.")
                    
                    # Normalize to desired size range
                    valid = size_col.notna()
                    sizes = np.full(len(hits_df), point_size_min)
                    
                    if valid.any():
                        values = size_col[valid].values
                        # Linear scaling from min to max value → min to max size
                        v_min, v_max = values.min(), values.max()
                        if v_max > v_min:
                            normalized = (values - v_min) / (v_max - v_min)
                        else:
                            normalized = np.ones_like(values) * 0.5
                        scaled = point_size_min + normalized * (point_size_max - point_size_min)
                        scaled *= point_size_scale
                        sizes[valid] = np.clip(scaled, point_size_min, point_size_max)
                    
                    # Assign sizes back
                    hits_df['_point_size'] = sizes
                    size_list = hits_df['_point_size'].tolist()
                    
                    print(f"[*] Hit point sizes scaled by '{point_size}' "
                        f"(range: {size_col.min():.2f} → {size_col.max():.2f})")
                else:
                    size_list = 40  # fixed size
                    if point_size:
                        print(f"[*] Warning: column '{point_size}' not found. Using fixed hit size.")
                
                # Plot hits with variable or fixed size
                sns.scatterplot(
                    data=hits_df, x='umap_x', y='umap_y',
                    color='red', alpha=0.9, s=size_list,
                    ax=ax, label='Hits', edgecolor='black', linewidth=0.5,
                    rasterized=True
                )
                ax.legend(markerscale=1.2, loc='upper right')
            
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.set_title(plot_title or 'Chemical Space UMAP')
            plt.tight_layout()
            pdf.savefig(fig, dpi=300)
            plt.close(fig)
        
        print(f"[*] Plot saved: {output_file}")
    
    def draw_hits(
        self,
        work_df: pd.DataFrame,
        max_hits: int = 30,
        mols_per_row: int = 5,
        sort_by: str = 'mw',
        ascending: bool = False,
        output_file: str = 'hit_structures.pdf',
        sub_img_size: tuple = (300, 300)
    ):
        hits = work_df[work_df['is_hit']].copy()
        if len(hits) == 0:
            return
        
        hits = hits.sort_values(sort_by, ascending=ascending).head(max_hits)
        mols = hits['mol'].tolist()
        
        # Prepare legends
        exclude_cols = [x for x in self.lib.df.columns if x not in ['mw', 'cluster']] + ['is_hit']
        legend_cols = [col for col in hits.columns if col not in exclude_cols]
        legends = []
        for _, row in hits.iterrows():
            lines = []
            for col in legend_cols:
                val = row[col]
                # Format special cases
                if pd.isna(val):
                    continue
                if isinstance(val, float):
                    if col.lower().startswith(('sim', 'tanimoto', 'similarity', 'score', 'prob')):
                        lines.append(f"{col}: {val:.3f}")
                    elif col.lower() in ('mw', 'molecular_weight', 'logp', 'tpsa'):
                        lines.append(f"{col}: {val:.1f}")
                    else:
                        lines.append(f"{col}: {val:.2f}")
                else:
                    # For strings, IDs, clusters, etc.
                    lines.append(f"{col}: {val}")
            
            # Always show compound_id and cluster near the top
            legend_text = f"ID: {row['compound_id']}"
            if lines:
                legend_text += "\n" + "\n".join(lines[:12])  # limit to ~12 lines to avoid overflow
            
            legends.append(legend_text)
        
        # for _, row in hits.iterrows():
        #     lines = [f"ID: {row['compound_id']}", f"Cluster: {row['cluster']}"]
        #     if 'tanimoto_similarity' in row:
        #         lines.append(f"Sim: {row['tanimoto_similarity']:.3f}")
        #     legends.append("\n".join(lines))
                
        with PdfPages(output_file) as pdf:
            i = 0
            while i < len(mols):
                page_mols = mols[i:i + mols_per_row]
                page_legends = legends[i:i + mols_per_row]
                n = len(page_mols)
                
                fig, axs = plt.subplots(1, mols_per_row, figsize=(mols_per_row*3, 4))
                if mols_per_row == 1: axs = [axs]
                
                for j in range(mols_per_row):
                    if j < n:
                        drawer = rdMolDraw2D.MolDraw2DCairo(sub_img_size[0], sub_img_size[1])
                        rdMolDraw2D.PrepareAndDrawMolecule(drawer, page_mols[j])
                        drawer.FinishDrawing()
                        img = mpimg.imread(BytesIO(drawer.GetDrawingText()), format='png')
                        axs[j].imshow(img)
                        axs[j].set_title(page_legends[j], fontsize=8)
                    axs[j].axis('off')
                
                pdf.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                i += mols_per_row
        print(f"[*] Structures saved: {output_file}")
        
    def export_results(self, work_df, output_file: str) -> None:
        os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else ".", exist_ok=True)
        export_cols = [x for x in work_df.columns if x not in ['fp','mol']]
        work_df[export_cols].to_excel(output_file, index=False)
        print(f"[*] Results exported to {output_file}")




# ==========================================
# MAIN EXECUTION BLOCK
# ==========================================
if __name__ == "__main__":

    # Configuration
    BASE_DIR = '/mnt/f/workspace/Project_p53'
    LIB_FILE = 'ASMS_primary/HTS-V5.0-370k.csv'
    N_NEIGHBORS = 10
    MIN_DIST = 0.05
    tody = datetime.now().strftime("%Y-%m-%d")
    # PKL_NAME = f'{tody}_library_n{N_NEIGHBORS}_d{MIN_DIST}.pkl'
    PKL_NAME = f'2026-01-10_library_n{N_NEIGHBORS}_d{MIN_DIST}.pkl'

    if not os.path.exists(BASE_DIR):
        print(f"Directory {BASE_DIR} not found. Check your paths.")
    else:
        os.chdir(BASE_DIR)

        # Step 1: Initialize or Load Library
        if os.path.exists(PKL_NAME):
            print(f"[*] Loading existing library from {PKL_NAME}")
            lib = joblib.load(PKL_NAME)
        else:
            lib = ChemicalLibraryAnalyzer()
            lib.load_library(LIB_FILE, id_col='Product No', smiles_col='smiles')
            lib.generate_fingerprints_and_properties()
            lib.select_diverse_representatives(num_reps=500)
            lib.run_umap(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST)
            joblib.dump(lib, PKL_NAME)
        
        # Step 2: Visualization
        viz = HitVisualizer(lib)
        
        # Case: Generic Library View
        if not os.path.exists(f'{tody}_library.xlsx'):
            df_lib = viz.add_hits(source=None)
            viz.map_umap(df_lib, output_file=f'{tody}_library.pdf')
            viz.export_results(df_lib, output_file=f'{tody}_library.xlsx')
        
        # Case: Specific Hits
        # hits_file = 'ASMS_primary/WuXi Lin Gang Lab ASMS-P53-0109_simple.csv'
        # sort_by = 'Height'
        hits_file = 'ASMS_confirm/WuXi Lin Gang Lab ASMS-P53-20260122_simple.csv'
        sort_by = 'RBA'
        if os.path.exists(hits_file) and (not os.path.exists(f'{tody}_hits_structure_{sort_by}.pdf')):
            df_hits = viz.add_hits(source=hits_file, id_col='Name')
            viz.map_umap(df_hits, output_file=f'{tody}_hits_umap_{sort_by}.pdf', point_size=sort_by)
            viz.export_results(df_hits, f'{tody}_hits_result_{sort_by}.xlsx')
            viz.draw_hits(df_hits, max_hits=300, sort_by='RBA', ascending=False, output_file=f'{tody}_hits_structure_{sort_by}.pdf')

        # # Case: Similarity Search
        # query_smiles = {'pc14586'="CNC(=O)C1=CC(=C(C=C1)NCC#CC2=CC3=C(C=CC=C3N2CC(F)(F)F)N[C@@H]4CCN(C[C@@H]4F)C)OC",
        #                 'activator7'="CNC(=O)C1=CC(=C(C=C1)NCC#CC2=CC3=C(C=CC=C3N2CC(F)(F)F)N[C@@H]4CCN(C[C@@H]4F)C)OC",
        #                 }
        # for k, v in query_smiles:
        #     if not os.path.exists(f'{tody}_similar_structure_{k}.pdf'):
        #         top_analogs = lib.search_similar(v, top_n=50)
        #         top_analogs = top_analogs[['compound_id','tanimoto_similarity']]
        #         df_search = viz.add_hits(source=top_analogs, id_col='compound_id', name='Analogs')
        #         viz.map_umap(df_search, output_file=f'{tody}_similar_umap_{k}.pdf', point_size='tanimoto_similarity')
        #         viz.export_results(df_search, f'{tody}_similar_result_{k}.xlsx')
        #         viz.draw_hits(df_search, max_hits=50, sort_by='tanimoto_similarity', ascending=False, output_file=f'{tody}_similar_structure_{k}.pdf')
       


