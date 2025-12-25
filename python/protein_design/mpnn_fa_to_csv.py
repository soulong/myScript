import os
import pandas as pd
import re
import argparse

def parse_fa_file(file_path):
    data = {'id': [], 'score': [], 'seq_recovery': [], 'sequence': []}
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    current_sub = None
    sample_counter = 0

    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('>sub_'):
            match = re.search(r'>sub_(\d+)', line)
            if match:
                current_sub = match.group(1)
                sample_counter = 0
        elif line.startswith('>T='):
            sample_counter += 1
            score_match = re.search(r'score=([\d.]+)', line)
            recovery_match = re.search(r'seq_recovery=([\d.]+)', line)
            score = float(score_match.group(1)) if score_match else None
            seq_recovery = float(recovery_match.group(1)) if recovery_match else None
            sequence_line = lines[i + 1].strip() if i + 1 < len(lines) else ''
            data['id'].append(f"sub_{current_sub}#{sample_counter}")
            data['score'].append(score)
            data['seq_recovery'].append(seq_recovery)
            data['sequence'].append(sequence_line)

    return pd.DataFrame(data)

def merge_fa_files(dir_path):
    all_data = []
    for root, _, files in os.walk(dir_path):
        for file in files:
            if file.endswith('.fa'):
                file_path = os.path.join(root, file)
                df = parse_fa_file(file_path)
                all_data.append(df)
    if all_data:
        return pd.concat(all_data, ignore_index=True)
    else:
        return pd.DataFrame(columns=['id', 'score', 'seq_recovery', 'sequence'])

def save_fasta(df, fasta_path):
    with open(fasta_path, 'w') as f:
        for i, row in df.iterrows():
            f.write(f">{row['id']}\n{row['sequence']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge .fa files to CSV and FASTA.")
    parser.add_argument("input_dir", help="Directory containing .fa files")
    parser.add_argument("output_csv", help="Output CSV file path")
    parser.add_argument("--output_fasta", help="Output FASTA file path")
    args = parser.parse_args()

    merged_df = merge_fa_files(args.input_dir)

    if not merged_df.empty:
        merged_df.to_csv(args.output_csv, index=False)
        print(f"[✓] CSV saved to: {args.output_csv}")
        if args.output_fasta is not None:
            save_fasta(merged_df, args.output_fasta)
            print(f"[✓] FASTA saved to: {args.output_fasta}")
    else:
        print("[!] No .fa files found in the specified directory.")
