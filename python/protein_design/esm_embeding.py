#!/usr/bin/env python3
"""
conda env: esm
Protein Sequence Embedding and Clustering Tool

This script/package uses ESM (Evolutionary Scale Modeling) to generate embeddings for protein sequences,
performs clustering using UMAP, and visualizes the results.

Usage as command line tool:
    python protein_clustering.py --input sequences.fasta --output results
    python protein_clustering.py --input sequences.csv --output results --format csv
    python protein_clustering.py --input sequences.fasta --model facebook/esm2_t12_35M_UR50D --clusters 5

Usage as Python package:
    from protein_clustering import ProteinEmbedder, ProteinClusterer, cluster_from_file
    
Requirements:
    pip install torch transformers umap-learn scikit-learn matplotlib seaborn pandas numpy h5py argparse
"""

import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from transformers import EsmModel, EsmTokenizer
import umap
from sklearn.cluster import KMeans, DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import h5py
import pickle
import argparse
import sys
import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Union
import warnings
warnings.filterwarnings('ignore')


class ProteinEmbedder:
    """Class to handle protein sequence embedding using ESM models"""
    
    def __init__(self, model_name: str = "facebook/esm2_t6_8M_UR50D", verbose: bool = True):
        """
        Initialize the protein embedder
        
        Args:
            model_name: ESM model name. Options include:
                - facebook/esm2_t6_8M_UR50D (smallest, fastest)
                - facebook/esm2_t12_35M_UR50D (medium)
                - facebook/esm2_t30_150M_UR50D (larger)
                - facebook/esm2_t33_650M_UR50D (large)
            verbose: Whether to print progress messages
        """
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.verbose = verbose
        
        if self.verbose:
            print(f"Using device: {self.device}")
            print(f"Loading ESM model: {model_name}...")
        
        self.tokenizer = EsmTokenizer.from_pretrained(model_name)
        self.model = EsmModel.from_pretrained(model_name)
        self.model.to(self.device)
        self.model.eval()
        
        if self.verbose:
            print("Model loaded successfully!")
    
    def embed_sequences(self, sequences: List[str], batch_size: int = 8) -> np.ndarray:
        """
        Generate embeddings for protein sequences
        
        Args:
            sequences: List of protein sequences
            batch_size: Batch size for processing
            
        Returns:
            Numpy array of embeddings (n_sequences, embedding_dim)
        """
        embeddings = []
        
        if self.verbose:
            print(f"Generating embeddings for {len(sequences)} sequences...")
        
        with torch.no_grad():
            for i in range(0, len(sequences), batch_size):
                batch_sequences = sequences[i:i+batch_size]
                if self.verbose:
                    print(f"Processing batch {i//batch_size + 1}/{(len(sequences)-1)//batch_size + 1}")
                
                # Tokenize sequences
                tokens = self.tokenizer(
                    batch_sequences, 
                    return_tensors="pt", 
                    padding=True, 
                    truncation=True,
                    max_length=1024  # ESM models have max length limits
                )
                
                # Move to device
                tokens = {k: v.to(self.device) for k, v in tokens.items()}
                
                # Get embeddings
                outputs = self.model(**tokens)
                
                # Use mean pooling over sequence length (excluding special tokens)
                attention_mask = tokens['attention_mask']
                token_embeddings = outputs.last_hidden_state
                
                # Mean pooling
                input_mask_expanded = attention_mask.unsqueeze(-1).expand(token_embeddings.size()).float()
                sum_embeddings = torch.sum(token_embeddings * input_mask_expanded, 1)
                sum_mask = torch.clamp(input_mask_expanded.sum(1), min=1e-9)
                batch_embeddings = sum_embeddings / sum_mask
                
                embeddings.append(batch_embeddings.cpu().numpy())
        
        return np.vstack(embeddings)


class ProteinClusterer:
    """Class to handle clustering and visualization of protein embeddings"""
    
    def __init__(self, verbose: bool = True):
        self.umap_model = None
        self.cluster_model = None
        self.embeddings = None
        self.umap_embeddings = None
        self.cluster_labels = None
        self.verbose = verbose
        
    def fit_umap(self, embeddings: np.ndarray, n_neighbors: int = 15, 
                 min_dist: float = 0.1, n_components: int = 2, 
                 metric: str = 'cosine') -> np.ndarray:
        """
        Fit UMAP dimensionality reduction
        
        Args:
            embeddings: High-dimensional embeddings
            n_neighbors: UMAP n_neighbors parameter
            min_dist: UMAP min_dist parameter
            n_components: Number of components for output
            metric: Distance metric
            
        Returns:
            Low-dimensional UMAP embeddings
        """
        if self.verbose:
            print("Fitting UMAP...")
        
        # Adjust n_neighbors if too large for dataset
        n_neighbors = min(n_neighbors, len(embeddings) - 1)
        
        self.umap_model = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            metric=metric,
            random_state=42
        )
        
        self.embeddings = embeddings
        self.umap_embeddings = self.umap_model.fit_transform(embeddings)
        
        if self.verbose:
            print(f"UMAP completed. Shape: {self.umap_embeddings.shape}")
        return self.umap_embeddings
    
    def cluster_sequences(self, method: str = 'kmeans', n_clusters: int = None, 
                         eps: float = 0.5, min_samples: int = 5) -> np.ndarray:
        """
        Cluster sequences based on UMAP embeddings
        
        Args:
            method: Clustering method ('kmeans' or 'dbscan')
            n_clusters: Number of clusters for K-means
            eps: DBSCAN eps parameter
            min_samples: DBSCAN min_samples parameter
            
        Returns:
            Cluster labels
        """
        if self.umap_embeddings is None:
            raise ValueError("Must fit UMAP first!")
        
        if self.verbose:
            print(f"Clustering using {method}...")
        
        if method == 'kmeans':
            if n_clusters is None:
                # Determine optimal number of clusters using silhouette score
                n_clusters = self._find_optimal_clusters()
            
            self.cluster_model = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            self.cluster_labels = self.cluster_model.fit_predict(self.umap_embeddings)
            
        elif method == 'dbscan':
            self.cluster_model = DBSCAN(eps=eps, min_samples=min_samples)
            self.cluster_labels = self.cluster_model.fit_predict(self.umap_embeddings)
            
        else:
            raise ValueError("Method must be 'kmeans' or 'dbscan'")
        
        n_clusters_found = len(np.unique(self.cluster_labels))
        if self.verbose:
            print(f"Found {n_clusters_found} clusters")
        
        return self.cluster_labels
    
    def _find_optimal_clusters(self, max_k: int = 10) -> int:
        """Find optimal number of clusters using silhouette score"""
        if len(self.umap_embeddings) < max_k:
            max_k = len(self.umap_embeddings) - 1
        
        if max_k < 2:
            return 2
        
        silhouette_scores = []
        k_range = range(2, min(max_k + 1, len(self.umap_embeddings)))
        
        for k in k_range:
            if k >= len(self.umap_embeddings):
                break
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = kmeans.fit_predict(self.umap_embeddings)
            score = silhouette_score(self.umap_embeddings, labels)
            silhouette_scores.append(score)
        
        if not silhouette_scores:
            return 2
        
        optimal_k = k_range[np.argmax(silhouette_scores)]
        if self.verbose:
            print(f"Optimal number of clusters: {optimal_k}")
        return optimal_k
    
    def plot_results(self, sequences: List[str] = None, sequence_names: List[str] = None,
                    figsize: Tuple[int, int] = (12, 8), save_path: str = None,
                    title: str = "Protein Sequence Clustering"):
        """
        Plot UMAP results with cluster colors
        
        Args:
            sequences: Original sequences (for hover info)
            sequence_names: Names/IDs for sequences
            figsize: Figure size
            save_path: Path to save the plot
            title: Plot title
        """
        if self.umap_embeddings is None or self.cluster_labels is None:
            raise ValueError("Must fit UMAP and cluster first!")
        
        plt.figure(figsize=figsize)
        
        # Create scatter plot
        unique_labels = np.unique(self.cluster_labels)
        colors = plt.cm.tab20(np.linspace(0, 1, len(unique_labels)))
        
        for i, label in enumerate(unique_labels):
            mask = self.cluster_labels == label
            plt.scatter(
                self.umap_embeddings[mask, 0], 
                self.umap_embeddings[mask, 1], 
                c=[colors[i]], 
                label=f'Cluster {label}',
                alpha=0.7,
                s=50
            )
        
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.title(title)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add cluster centers for K-means
        if hasattr(self.cluster_model, 'cluster_centers_'):
            centers_umap = self.cluster_model.cluster_centers_
            plt.scatter(centers_umap[:, 0], centers_umap[:, 1], 
                       c='red', marker='x', s=200, linewidths=3, 
                       label='Centroids', alpha=0.8)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            if self.verbose:
                print(f"Plot saved to {save_path}")
        
        return plt.gcf()


def load_sequences_from_fasta(fasta_path: str) -> Tuple[List[str], List[str]]:
    """
    Load sequences from FASTA file
    
    Args:
        fasta_path: Path to FASTA file
        
    Returns:
        Tuple of (sequences, sequence_names)
    """
    sequences = []
    names = []
    current_seq = ""
    current_name = ""
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                    names.append(current_name)
                current_name = line[1:]  # Remove '>'
                current_seq = ""
            else:
                current_seq += line
        
        # Add last sequence
        if current_seq:
            sequences.append(current_seq)
            names.append(current_name)
    
    return sequences, names


def load_sequences_from_csv(csv_path: str, id_col: str = None, seq_col: str = None) -> Tuple[List[str], List[str]]:
    """
    Load sequences from CSV file
    
    Args:
        csv_path: Path to CSV file
        id_col: Name of ID column (default: first column)
        seq_col: Name of sequence column (default: second column)
        
    Returns:
        Tuple of (sequences, sequence_names)
    """
    df = pd.read_csv(csv_path)
    
    # Use column names or indices
    if id_col is None:
        id_col = df.columns[0]
    if seq_col is None:
        seq_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
    
    sequences = df[seq_col].astype(str).tolist()
    names = df[id_col].astype(str).tolist()
    
    return sequences, names


def save_results_h5(embeddings: np.ndarray, umap_embeddings: np.ndarray, 
                   cluster_labels: np.ndarray, sequences: List[str], 
                   sequence_names: List[str], output_path: str):
    """
    Save all results to HDF5 file
    
    Args:
        embeddings: Original high-dimensional embeddings
        umap_embeddings: UMAP embeddings
        cluster_labels: Cluster labels
        sequences: Original sequences
        sequence_names: Sequence names/IDs
        output_path: Output HDF5 file path
    """
    with h5py.File(output_path, 'w') as f:
        # Save embeddings
        f.create_dataset('embeddings', data=embeddings)
        f.create_dataset('umap_embeddings', data=umap_embeddings)
        f.create_dataset('cluster_labels', data=cluster_labels)
        
        # Save sequences and names as string arrays
        seq_dt = h5py.special_dtype(vlen=str)
        f.create_dataset('sequences', data=sequences, dtype=seq_dt)
        f.create_dataset('sequence_names', data=sequence_names, dtype=seq_dt)
        
        # Save metadata
        f.attrs['n_sequences'] = len(sequences)
        f.attrs['embedding_dim'] = embeddings.shape[1]
        f.attrs['n_clusters'] = len(np.unique(cluster_labels))


def save_results_csv(embeddings: np.ndarray, umap_embeddings: np.ndarray, 
                    cluster_labels: np.ndarray, sequences: List[str], 
                    sequence_names: List[str], output_path: str):
    """
    Save results to CSV file
    
    Args:
        embeddings: Original high-dimensional embeddings
        umap_embeddings: UMAP embeddings
        cluster_labels: Cluster labels
        sequences: Original sequences
        sequence_names: Sequence names/IDs
        output_path: Output CSV file path
    """
    df_data = {
        'sequence_name': sequence_names,
        'sequence': sequences,
        'cluster': cluster_labels,
        'umap_1': umap_embeddings[:, 0],
        'umap_2': umap_embeddings[:, 1]
    }
    
    # Add embedding dimensions
    for i in range(embeddings.shape[1]):
        df_data[f'embedding_{i}'] = embeddings[:, i]
    
    df = pd.DataFrame(df_data)
    df.to_csv(output_path, index=False)


def cluster_from_file(input_path: str, output_prefix: str = "results", 
                     input_format: str = "auto", model_name: str = "facebook/esm2_t6_8M_UR50D",
                     n_clusters: int = None, method: str = "kmeans", 
                     batch_size: int = 8, verbose: bool = True) -> Dict:
    """
    Complete clustering pipeline from file input
    
    Args:
        input_path: Path to input file (FASTA or CSV)
        output_prefix: Prefix for output files
        input_format: Input format ('fasta', 'csv', or 'auto')
        model_name: ESM model name
        n_clusters: Number of clusters (auto-detect if None)
        method: Clustering method ('kmeans' or 'dbscan')
        batch_size: Batch size for embedding generation
        verbose: Whether to print progress
        
    Returns:
        Dictionary with results
    """
    
    # Determine input format
    if input_format == "auto":
        if input_path.lower().endswith(('.fasta', '.fa', '.fas')):
            input_format = "fasta"
        elif input_path.lower().endswith('.csv'):
            input_format = "csv"
        else:
            raise ValueError("Cannot determine input format. Please specify 'fasta' or 'csv'")
    
    # Load sequences
    if verbose:
        print(f"Loading sequences from {input_path} (format: {input_format})")
    
    if input_format == "fasta":
        sequences, names = load_sequences_from_fasta(input_path)
    elif input_format == "csv":
        sequences, names = load_sequences_from_csv(input_path)
    else:
        raise ValueError("Input format must be 'fasta' or 'csv'")
    
    if verbose:
        print(f"Loaded {len(sequences)} sequences")
    
    # Generate embeddings
    embedder = ProteinEmbedder(model_name=model_name, verbose=verbose)
    embeddings = embedder.embed_sequences(sequences, batch_size=batch_size)
    
    # Clustering
    clusterer = ProteinClusterer(verbose=verbose)
    umap_embeddings = clusterer.fit_umap(embeddings)
    cluster_labels = clusterer.cluster_sequences(method=method, n_clusters=n_clusters)
    
    # Save results
    if verbose:
        print("Saving results...")
    
    # Save HDF5
    h5_path = f"{output_prefix}.h5"
    save_results_h5(embeddings, umap_embeddings, cluster_labels, sequences, names, h5_path)
    
    # Save CSV
    csv_path = f"{output_prefix}.csv"
    save_results_csv(embeddings, umap_embeddings, cluster_labels, sequences, names, csv_path)
    
    # Save plot as PDF
    pdf_path = f"{output_prefix}_plot.pdf"
    fig = clusterer.plot_results(sequences, names, title=f"Protein Clustering - {Path(input_path).name}")
    
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    
    if verbose:
        print(f"Results saved:")
        print(f"  - {h5_path} (HDF5 format)")
        print(f"  - {csv_path} (CSV format)")
        print(f"  - {pdf_path} (Plot)")
    
    plt.close(fig)
    
    # Return results
    return {
        'embeddings': embeddings,
        'umap_embeddings': umap_embeddings,
        'cluster_labels': cluster_labels,
        'sequences': sequences,
        'sequence_names': names,
        'n_clusters': len(np.unique(cluster_labels))
    }



def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description="Protein Sequence Embedding and Clustering Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with FASTA file
  python protein_clustering.py --input sequences.fasta --output results
  
  # CSV input with custom parameters
  python protein_clustering.py --input sequences.csv --format csv --output results --clusters 5
  
  # Using larger model
  python protein_clustering.py --input sequences.fasta --model facebook/esm2_t12_35M_UR50D --output results
  
  # DBSCAN clustering
  python protein_clustering.py --input sequences.fasta --method dbscan --output results
        """
    )
    
    parser.add_argument('--input', '-i', required=True, help='Input file (FASTA or CSV)')
    parser.add_argument('--output', '-o', default='results', help='Output prefix (default: results)')
    parser.add_argument('--format', '-f', choices=['fasta', 'csv', 'auto'], default='auto', 
                       help='Input format (default: auto-detect)')
    parser.add_argument('--model', '-m', default='facebook/esm2_t6_8M_UR50D',
                       help='ESM model name (default: facebook/esm2_t6_8M_UR50D)')
    parser.add_argument('--clusters', '-c', type=int, default=None,
                       help='Number of clusters (default: auto-detect)')
    parser.add_argument('--method', choices=['kmeans', 'dbscan'], default='kmeans',
                       help='Clustering method (default: kmeans)')
    parser.add_argument('--batch-size', '-b', type=int, default=8,
                       help='Batch size for embedding generation (default: 8)')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Suppress progress messages')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)
    
    try:
        # Run clustering pipeline
        results = cluster_from_file(
            input_path=args.input,
            output_prefix=args.output,
            input_format=args.format,
            model_name=args.model,
            n_clusters=args.clusters,
            method=args.method,
            batch_size=args.batch_size,
            verbose=not args.quiet
        )
        
        # Print summary
        if not args.quiet:
            print(f"\n=== Analysis Summary ===")
            print(f"Input file: {args.input}")
            print(f"Number of sequences: {len(results['sequences'])}")
            print(f"Embedding dimension: {results['embeddings'].shape[1]}")
            print(f"Number of clusters: {results['n_clusters']}")
            print(f"Cluster distribution: {np.bincount(results['cluster_labels'])}")
            print(f"\nAnalysis completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        if not args.quiet:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()