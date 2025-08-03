
# this is used for convert protein designed

import pandas as pd
import sys

def table_to_fasta(table_file, output_fasta, seq_colname="description", fasta_colname="sequence"):
  """Converts a table file with sequence names and sequences to a FASTA file.

  Args:
    table_file: The path to the input table file.
    output_fasta: The path to the output FASTA file.
    seq_colname: colname of sequence name
    fasta_colname: colname of AA sequence
  """

  # Read the table file into a pandas DataFrame
  df = pd.read_csv(table_file)

  # Create a list to store the FASTA sequences
  fasta_sequences = []

  # Iterate over the rows of the DataFrame
  for index, row in df.iterrows():
    seq_name = row[seq_colname]  # Replace 'seq_name' with the actual column name
    protein_sequence = row[fasta_colname]  # Replace 'protein_sequence' with the actual column name
    fasta_sequence = f">{seq_name}\n{protein_sequence}"
    fasta_sequences.append(fasta_sequence)

  # Write the FASTA sequences to a file
  with open(output_fasta, 'w') as f:
    f.write('\n'.join(fasta_sequences))


if __name__ == "__main__":
  if len(sys.argv) != 3:
    print("Usage: python table_to_fasta.py <input_table_file> <output_fasta_file>")
    sys.exit(1)

  table_file = sys.argv[1]
  output_fasta = sys.argv[2]
  table_to_fasta(table_file, output_fasta)