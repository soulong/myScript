import os
import requests
import re
import shutil
import pandas as pd

from pathlib import Path
from tqdm import tqdm



def msa_search(sequence, API_KEY, databases=['Uniref30_2302', 'colabfold_envdb_202108', 'PDB70_220313']):
    msa_search_url = "https://health.api.nvidia.com/v1/biology/colabfold/msa-search/predict"
    payload = {
        "sequence": sequence,
        "databases": databases,
        "e_value": 0.0001,
        "iterations": 1,
        "max_msa_sequences": 10000,
        "run_structural_template_search": False,
        "output_alignment_formats": ["a3m"],
    }
    headers = {
        "Authorization": f"Bearer {API_KEY}",
        "content-type": "application/json",
        "NVCF-POLL-SECONDS": "300",
    }
    # Call MSA-Search NIM
    response = requests.post(msa_search_url, json=payload, headers=headers)
    return response.json()


def parse_sequences(input_string, n, sequence):
    """
    Parse the output of alignments from the MSA-Search NIM to be used downstream

    Args:
        input_string (str): The output file of alignments in a string format
        n (int): The amount of alignments to return from the output when parsing
        sequence (str): The query sequence for alignment

    Returns:
        list: A list of alignment identifiers and sequences, starting with the query,
              where the amount of sequences is given by n
    """
    # Output is parsed to have a line for the sequence id and sequence itself so `n` returns correlates to n*2 lines
    n = n * 2
    # First, handle the `Query` block separately
    lines = input_string.strip().split('\n')
    # Now process the rest of the lines
    remaining_string = "\n".join(lines[:])
    # Regex to find blocks starting with `>` and then followed by a sequence.
    pattern = re.compile(r'\n>(.*?)\n(.*?)(?=\n>|\Z)', re.DOTALL)
    matches = pattern.finditer(remaining_string)
    output_list_to_order = []
    for match in matches:
        # The name is the first capturing group, split by tab and take the first part
        name_full = match.group(1).split('\t')[0]
        SW_score = match.group(1).split('\t')[1]
        # The sequence is the second capturing group
        sequence_raw = match.group(2).strip()
        aligned_sequence = ''.join(char for char in sequence_raw if char.isupper() or not char.isalpha())
        # Store the aligned sequence in the list of outputs by name, sequence, Smith-Waterman score
        output_list_to_order.append((f'>{name_full}', aligned_sequence, int(SW_score)))
    output_lines = output_list_to_order[:n]
    return output_lines


def validate_a3m_format(alignments_string):
    """
    Validate that the alignment string follows A3M format.

    Args:
        alignments_string (str): String containing alignments

    Returns:
        bool: True if valid A3M format, False otherwise
    """
    lines = alignments_string.strip().split('\n')
    if len(lines) < 2:
        return False

    # Check that we have alternating header and sequence lines
    for i, line in enumerate(lines):
        if i % 2 == 0:  # Even indices should be headers
            if not line.startswith('>'):
                return False
        else:  # Odd indices should be sequences
            if line.startswith('>'):
                return False
            # Sequences should only contain valid amino acid characters and gaps
            if not all(c in 'ACDEFGHIKLMNPQRSTVWY-' for c in line.upper()):
                return False

    return True


def write_alignments_to_a3m(alignments_data, uniprot_id, output_dir):
    """
    Write alignment data to a3M format file.

    Args:
        alignments_data: Either a list of alternating headers/sequences or a string containing alignments
        uniprot_id (str): Uniprot ID of the protein
        output_dir (str): Directory for the output a3M file

    Returns:
        str: Path to the created a3M file
    """
    # Ensure output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    output_path = Path(output_dir) / f"{uniprot_id}_msa_alignments.a3m"

    # Handle both list and string input formats
    if isinstance(alignments_data, list):
        alignments_string = '\n'.join(alignments_data)
    elif isinstance(alignments_data, str):
        alignments_string = alignments_data
    else:
        raise ValueError("alignments_data must be either a list or string")

    # Validate A3M format
    if not alignments_string.strip():
        raise ValueError("Empty alignment data provided")

    # Count sequences for reporting
    sequence_count = alignments_string.count('>')
    if sequence_count == 0:
        raise ValueError("No sequences found in alignment data")

    # Validate A3M format structure
    if not validate_a3m_format(alignments_string):
        print("Warning: Alignment data may not follow strict A3M format")
        print("Proceeding with file creation...")

    print(f"Writing {sequence_count} sequences to A3M format: {output_path}")

    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            # Write the alignments
            f.write(alignments_string)
            # Ensure file ends with newline
            if not alignments_string.endswith('\n'):
                f.write('\n')

        # Verify the file was created successfully
        if output_path.exists():
            file_size = output_path.stat().st_size
            print(f"Successfully created A3M file:")
            print(f"File: {output_path}")
            print(f"Size: {file_size:,} bytes")
            print(f"Sequences: {sequence_count}")

            # Download the file to the user's machine
            try:
                files.download(str(output_path))
                print(f"File downloaded successfully: {output_path}")
            except Exception as download_error:
                print(f"Warning: Could not download file automatically: {download_error}")
                print(f"File is available at: {output_path}")

            return str(output_path)
        else:
            raise IOError(f"Failed to create file {output_path}")
    except Exception as e:
        print(f"Error writing A3M file: {e}")
        raise



def process_msa_alignments(msa_response_dict, sequence, uniprot_id, output_dir, databases=['Uniref30_2302', 'colabfold_envdb_202108', 'PDB70_220313'], max_sequences_per_db=10000):
    """
    Process MSA alignments from multiple databases and merge them into A3M format.

    Args:
        msa_response_dict (dict): MSA response data containing alignments
        sequence (str): Query sequence for alignment
        uniprot_id (str): Uniprot ID of the protein
        output_dir (str): Output directory for the A3M file
        databases (list): List of database names to process
        max_sequences_per_db (int): Maximum number of sequences to parse per database

    Returns:
        str: Path to the created A3M file
    """
    all_parsed_dataset_output = []
    for database in databases:
        print(f"Parsing results from database: {database}")
        # Pull string of alignments stored in json output for specific dataset
        a3m_dict_msa_search = msa_response_dict['alignments'][database]['a3m']['alignment']
        a3m_dict_msa_search_parsed = parse_sequences(a3m_dict_msa_search, max_sequences_per_db, sequence)
        num_sequences_aligned = (len(a3m_dict_msa_search_parsed))
        print(f"Number of sequences aligned: {num_sequences_aligned}")
        all_parsed_dataset_output.extend(a3m_dict_msa_search_parsed)
    # Sort all the alignments based off of the alignment score
    all_parsed_dataset_output.sort(key=lambda x: x[2], reverse=True)
    # Now that the alignments across all datasets are sorted, reformat each entry to name and sequence
    sorted_parsed_output_formatted = []
    for align_tuple in all_parsed_dataset_output:
        sorted_parsed_output_formatted.append(align_tuple[0])
        sorted_parsed_output_formatted.append(align_tuple[1])
    merged_alignments_protein = [f">query_sequence\n{sequence}"]
    merged_alignments_protein.extend(sorted_parsed_output_formatted)
    print(f"Total merged alignments: {len(merged_alignments_protein)}")
    # Write merged_alignments_protein to a3M format
    a3m_file_path = write_alignments_to_a3m(
        merged_alignments_protein,
        uniprot_id,
        output_dir
    )
    return a3m_file_path



if __name__ == 'main':

    API_KEY = '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'
    MSA_DATABASES = ['Uniref30_2302', 'colabfold_envdb_202108', 'PDB70_220313']

    output_dir = "./"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        os.makedirs(output_dir)

    csv_file_path = "target.csv"
    df = pd.read_csv(csv_file_path, low_memory=False)
    print(df.shape)
    df.head()

    sequences = tuple(zip(df['id'].tolist(), df['sequence'].tolist()))
    sequences = sorted(list(set(sequences)))
    print(f"Number of unique sequences: {len(sequences)} from parent dataset of {df.shape[0]} sequences")
    sequences[0]
    
    # Retrieve msa from nivida-nim
    for seq_id, seq in tqdm(sequences):
        try:
            print(f"\nProcessing protein: {seq_id}")
            print(f"Sequence length: {len(seq)}")

            # Call MSA-Search NIM
            msa_response_dict = msa_search(seq, API_KEY, MSA_DATABASES)

            # Check if the response contains the expected data
            if 'alignments' not in msa_response_dict:
                print(f"Warning: No alignments found for {seq_id}")
                continue

            # Process and create A3M file
            a3m_file_path = process_msa_alignments(msa_response_dict, seq, seq_id, output_dir)
            print(f"Successfully processed {seq_id} -> {a3m_file_path}")

        except Exception as e:
            print(f"Error processing {seq_id}: {e}")
            continue


    # List all created A3M files
    import glob
    a3m_files = glob.glob(f"{output_dir}/*.a3m")
    a3m_files = sorted(a3m_files)
    print(f"Created {len(a3m_files)} A3M files:")
    for file_path in a3m_files:
        file_size = Path(file_path).stat().st_size
        print(f"  - {Path(file_path).name} ({file_size:,} bytes)")

    print(f"\nAll A3M files are available in: {output_dir}")
    print("Files have been automatically downloaded to your machine.")