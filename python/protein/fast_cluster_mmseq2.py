import argparse
import subprocess

def cluster_fasta_with_mmseq2(fasta_file, output_prefix, min_seq_id, c, cov_mode):
    """
    Clusters protein sequences in a FASTA file using mmseq2.

    Args:
        fasta_file: Path to the FASTA file containing the protein sequences.
        output_prefix: Prefix for the output files.
        min_seq_id: Minimum sequence identity threshold for clustering.
        min_cov: Minimum sequence coverage threshold for clustering.
    """

    # Set mmseq2 command and options
    mmseq2_cmd = ["mmseqs", "easy-cluster", fasta_file, output_prefix, "tmp", "--min-seq-id", str(min_seq_id), "-c", str(c), "--cov-mode", str(cov_mode)]

    # Run mmseq2
    try:
        subprocess.run(mmseq2_cmd, check=True)
        print("Clustering completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running mmseq2: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster protein sequences using mmseq2.")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    parser.add_argument("output_prefix", help="Prefix for the output files")
    parser.add_argument("--min-seq-id", type=float, default=0.5, help="Minimum sequence identity threshold (default: 0.5)")
    parser.add_argument("-c", type=float, default=0.8, help="Minimum sequence coverage threshold (default: 0.8)")
    parser.add_argument("--cov-mode", type=int, default=1, help="Cover mode, int from 0 to 5")
    args = parser.parse_args()

    cluster_fasta_with_mmseq2(args.fasta_file, args.output_prefix, args.min_seq_id, args.c, args.cov_mode)