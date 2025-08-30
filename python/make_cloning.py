#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 18:22:45 2025

Make plasmid from vector and insert fasta files.

Should run in "ngs" conda environment

@author: hao
"""

#%%
import os, re, string, math
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch as Batch
import autosnapgene as snap
# from google_crc32c import exc
from tqdm import tqdm
from typing import List, AnyStr, Optional
import pandas as pd


def read_dna_file(f: Optional[AnyStr] = None, multiple: bool = False) -> Seq | List[Seq]:
    '''
    read .fa or .dna file, 
    usually a Bio.Seq.Seq object return, if .fa contain multuple records and multiple=True, will return a list of Seq object
    
    # example: read_dna('ZQ360_Lenti_G11-Xba1-shuttle-BamH1-NLS-Gal4DBD-HT3-T2A-G1-9.fa')
    '''
    if f.endswith('.fa'):
        data = list(SeqIO.parse(f, 'fasta'))
        if len(data) > 1 and multiple:
            print('multiple records found')
            seq = [x.seq for x in data]
        else:
            seq = data[0].seq

    elif f.endswith('.dna'):
        data = snap.parse(f)
        # feats = data.features
        seq = Seq(data.sequence)

    else:
        print('input must be .fa or .dna')
        seq = None

    return seq.upper()



def split_DNA_by_sequence(dna: Seq | AnyStr, split_sequence: List[Seq | AnyStr]) -> List[Seq]:
    '''
    split dna sequence by giving sequence list
    sequence must be find unique across whole sequence
    return Seq list, separated by split_sequence, length = len(split_sequence) + 1

    # split_DNA_by_sequence('AGGCTAACGCTTACCCGGGATTTTTGGGGG', ['GCTAA', 'CCGG', "TTTTT"])
    '''
    assert type(split_sequence) is list, 'split_sequence must be a list'
    split_sequence = [x.upper() for x in split_sequence]

    if type(dna) is str:
        dna = Seq(dna)
    dna = dna.upper()

    index = [index for index, _ in dna.search(split_sequence)]
    print(f'split index at {index}')
    assert len(index) == len(split_sequence), 'split_sequence is not unique in input dna sequence'
    
    splits = []
    prev = 0
    split_sequence_len = [len(x) for x in split_sequence]
    index_include_sequence = [x + y for x, y in zip(index, split_sequence_len)]
    for idx, num in enumerate(index):
        split = dna[prev:num]
        splits.append(split)
        prev = index_include_sequence[idx]
    splits.append(dna[prev:])
    
    return splits



def split_DNA_by_enzyme(dna: Seq | AnyStr, enzyme: List[AnyStr]) -> List[Seq]:
    '''
    split dna sequence by giving enzyme list
    sequence must be find unique across whole sequence
    return Seq list, this will generate segments seamlessly, length = len(enzyme) + 1

    # split_DNA_by_enzyme('CAGCGGTGGAGGCTCTAGAATGGAGGAGGCGGAGTATATGAATATGGCGCCGCAGGGATCCAGTG', ['XbaI', 'BamHI'])
    '''
    assert type(enzyme) is list, 'split_sequence must be a list'

    if type(dna) is str:
        dna = Seq(dna)
    dna = dna.upper()

    rb = Batch(enzyme)
    index = []
    # enzyme = ['XbaI', 'BamHI']
    for x in enzyme:
        res = rb.get(x).search(dna)
        assert len(res) == 1, f'{x} recognize site is not unique or not found'
        index.append(res[0] - 1)
    print(f'enzyme cut site at: {index}')

    splits = []
    prev = 0
    for _, num in enumerate(index):
        splits.append(dna[prev:num])
        prev = num
    splits.append(dna[prev:len(dna)])

    return splits




def make_plasmid(vector_file: AnyStr, # Path to the vector file, should be a .dna file
                insert_dir: Optional[AnyStr] = None, # Directory path containing insert files with .fa or .dna extensions
                insert_name_regex: Optional[AnyStr] = None, # Regex pattern to extract insert names from filenames, used when insert_dir is provided
                insert_sequences: Optional[List[AnyStr]] = None, # List of DNA sequences to be used as inserts
                insert_names: Optional[List[AnyStr]] = None, # List of names corresponding to insert_sequences, used when insert_sequences is provided
                insert_is_clean: bool = True, # Boolean to indicate if inserts need cleavage (removing sides by HR_sequence or enzyme)
                HR_sequence: Optional[List[AnyStr]] = None, # List of homology region sequences for splitting vector/insert, e.g., ['GATCTGGCAG', 'GTGGCGG']
                enzyme: Optional[List[AnyStr]] = None, # List of restriction enzyme names for splitting vector/insert, e.g., ['XbaI', 'BamHI']
                primer_overhang_vector_len: int = 5, # Length of vector overhang in primers, usually 25+ for HR, 4+ for enzyme-based
                primer_overhang_insert_len: int = 25 # Length of insert overhang in primers, used for PCR priming
                ) -> dict:
    """
    Constructs plasmids by combining a vector with insert sequences, either from a directory of files or provided sequences.
    The function processes the vector and inserts based on homology regions (HR) or restriction enzymes, generates constructs,
    and designs primers for PCR amplification.

    Args:
        vector_file (AnyStr): Path to the vector file, expected to be in .dna format.
        insert_dir (AnyStr, optional): Directory containing insert files (.fa or .dna). Defaults to None.
        insert_name_regex (AnyStr, optional): Regex pattern to extract names from insert filenames. Defaults to None.
        insert_sequences (List[AnyStr], optional): List of DNA sequences to use as inserts. Defaults to None.
        insert_names (List[AnyStr], optional): Names corresponding to insert_sequences. Defaults to None.
        insert_is_clean (bool, optional): If True, inserts are used as-is; if False, they are cleaved based on HR_sequence or enzyme. Defaults to True.
        HR_sequence (List[AnyStr], optional): Homology region sequences for splitting, e.g., ['GATCTGGCAG', 'GTGGCGG']. Defaults to None.
        enzyme (List[AnyStr], optional): Restriction enzyme names for splitting, e.g., ['XbaI', 'BamHI']. Defaults to None.
        primer_overhang_vector_len (int, optional): Length of vector overhang in primers, typically 25+ for HR. Defaults to 25.
        primer_overhang_insert_len (int, optional): Length of insert overhang in primers for PCR. Defaults to 25.

    Returns:
        dict: A dictionary containing:
            - constructs (List[Seq]): List of constructed plasmid sequences.
            - names (List[str]): List of names for the constructs.
            - insert_seqs (List[Seq]): List of processed insert sequences.
            - vector_seq (Seq): The original vector sequence.
            - vector_feats (List): Features of the vector.
            - primers (pd.DataFrame): DataFrame with primer information including forward, reverse, and reverse complement sequences.
    """

    assert os.path.exists(vector_file), "vector_file not existed"
    assert ((HR_sequence is None) + (enzyme is None )) == 1, 'must one HR_sequence or enzyme is not None'
    assert ((insert_dir is None) + (insert_sequences is None )) == 1, 'must one insert_dir or insert_sequences is not None'
    if HR_sequence:
        assert type(HR_sequence) is list, 'HR_sequence must be a list'
    if enzyme:
        assert type(enzyme) is list, 'enzyme must be a list'
    
    # vector read
    vector = snap.parse(vector_file)
    vector_feats = vector.features
    vector_seq = Seq(vector.sequence)

    # vector split
    print('------------- process vector -------------')
    process_func = split_DNA_by_sequence if HR_sequence else split_DNA_by_enzyme
    params = HR_sequence if HR_sequence else enzyme
    vector_splits = process_func(vector_seq, params)
    print('vector splits', *vector_splits, sep='\n')


    # read insert input
    print('------------- read inserts -------------')
    if insert_dir:
        print('use insert_dir as insertion')
        insert_files = glob(insert_dir + '/*.fa', recursive=True) + glob(insert_dir + '/*.dna', recursive=True)
        assert len(insert_files) > 0, 'no files found under insert_dir'
        insert = [read_dna_file(f, multiple=False) for f in insert_files]
        insert_names = [os.path.basename(f) for f in insert_files]
        print(f'total inserts [{len(insert_files)}] >>>', *insert_files, sep='\n')
    else:
        print('use insert_sequences as insertion')
        insert = [Seq(x) for x in insert_sequences]
        if insert_names:
            assert len(insert_sequences) == len(insert_names), 'insert_sequences and insert_names should have same length'
        else:
            insert_names = [x[:20] for x in insert_sequences] # use part of sequence as name
        print(f'total inserts [{len(insert_names)}] >>>', *insert_names, sep='\n')

    # process insert input
    if not insert_is_clean:
        print('------------- process inserts -------------')
        insert_processed = []
        insert_names_processed = []
        for n, x in zip(insert_names, insert):
            try:
                insert_use = process_func(x, params)[1] # use middle part as insert
                insert_processed.append(insert_use)
                insert_names_processed.append(n)
                print(f'insert length: {len(insert_use)}')
                print(f'{n}\n')
            except Exception as e:
                print(e)
        print(f'prcoessed inserts: {len(insert_names_processed)}')
    else:
        insert_processed = insert
        insert_names = insert_names
    


    print('------------- process constructs -------------')
    constructs = []
    names = []
    primers = pd.DataFrame(columns=['name','primer_f','primer_r','primer_r_revcomp'])
    for n, x in zip(insert_names, insert_processed):

        if insert_name_regex is not None:
            try:
                fname = re.search(insert_name_regex, n).group(1)
            except Exception as e:
                print(e)
                fname = n
        else:
            fname = n
        names.append(fname)

        vector_left = vector_splits[0] + HR_sequence[0] if HR_sequence else vector_splits[0]
        vector_right = HR_sequence[1] + vector_splits[2] if HR_sequence else vector_splits[2]
        constructs.append(vector_left + x + vector_right)

        # get primers
        primer_f = vector_left[-primer_overhang_vector_len:].lower() + x[:primer_overhang_insert_len]
        primer_r = x[-primer_overhang_insert_len:] + vector_right[:primer_overhang_vector_len].lower()
        primer_r_revcomp = primer_r.reverse_complement()
        primers.loc[len(primers)] = [fname, primer_f, primer_r, primer_r_revcomp]


    return {'constructs':constructs, 'names':names, 'insert_seqs':insert_processed,
            'vector_seq':vector_seq, 'vector_feats':vector_feats, 
            'primers':primers}




def write_snapgene(construct: Seq | AnyStr, 
                   path: AnyStr,
                   insert_seq: Seq | AnyStr, insert_name: AnyStr, insert_color: str = '#31849b',
                   vector_seq: Seq | AnyStr = None, vector_feats = None):
    construct = str(construct)
    
    tac = snap.SnapGene()
    tac.sequence = construct
    tac.set_topology('circular')

    # add vector features
    if vector_seq is not None and vector_feats is not None:
        for feat in vector_feats:
            # feat = vector_feats[2]
            try:
                update_seq = str(vector_seq[feat.segments[0].range[0]-1:feat.segments[0].range[1]])
                if construct.find(update_seq):
                    update_feat = snap.Feature.from_segment(
                        name=feat.name, 
                        color=feat.segments[0].color, 
                        # directionality = feat.directionality,
                        type=feat.type)
                    tac.add_feature(update_feat, update_seq)
            except:
                continue

    # add insert features   
    new_feat = snap.Feature.from_segment(name=insert_name, type='CDS', color=insert_color)
    tac.add_feature(new_feat, str(insert_seq))

    # save
    tac.write(path)

    return None



def convert_to_plate_format(primer, plate_nrow=8, plate_ncol=12):
    """
    Convert primer data to plate format, creating multiple plates if needed.
    
    Parameters:
    primer: DataFrame with primer data
    plate_nrow: number of rows per plate (default 8)
    plate_ncol: number of columns per plate (default 12)
    
    Returns:
    List of DataFrames, one for each plate needed
    """
    
    # Calculate plate capacity and number of plates needed
    plate_capacity = plate_nrow * plate_ncol
    num_primers = len(primer)
    num_plates = math.ceil(num_primers / plate_capacity)
    
    print(f"Number of primers: {num_primers}")
    print(f"Plate capacity: {plate_capacity}")
    print(f"Number of plates needed: {num_plates}")
    
    # List to store all plate DataFrames
    all_plates = []
    
    for plate_num in range(num_plates):
        # Calculate start and end indices for this plate
        start_idx = plate_num * plate_capacity
        end_idx = min((plate_num + 1) * plate_capacity, num_primers)
        
        # Subset primers for this plate
        primer_subset = primer.iloc[start_idx:end_idx].copy()
        primer_subset['id'] = range(1, len(primer_subset) + 1)
        
        # Create well names for this plate
        rows = list(string.ascii_uppercase[:plate_nrow]) * plate_ncol
        cols = [col for col in range(1, plate_ncol + 1) for _ in range(plate_nrow)]
        wells = [f"{row}{col}" for row, col in zip(rows, cols)]
        
        plate = pd.DataFrame({'well': wells})
        plate['id'] = range(1, len(plate) + 1)
        
        # Merge plate and primer data
        plate_primer = pd.merge(plate, primer_subset, on='id', how='outer')
        plate_primer = plate_primer.dropna(subset=['name'])
        
        # Convert to long shape
        plate_primer_processed = plate_primer.drop(columns=['id', 'primer_r'])
        plate_primer_processed = plate_primer_processed.rename(columns={
            'primer_f': '_f', 'primer_r_revcomp': '_r'})
        
        # Pivot longer
        plate_primer_long = pd.melt(
            plate_primer_processed,
            id_vars=['well', 'name'], value_vars=['_f', '_r'],
            var_name='order', value_name='sequence')
        
        # Unite name and order columns
        plate_primer_long['name'] = plate_primer_long['name'] + plate_primer_long['order']
        plate_primer_long = plate_primer_long.drop(columns=['order'])
        
        # Sort by natural order (A1, A2, ..., A9, A10, A11, A12) and name
        def natural_sort_key(well):
            """Create sort key for natural ordering of wells"""
            match = re.match(r'([A-Z])(\d+)', well)
            if match:
                letter, number = match.groups()
                return (letter, int(number))
            return (well, 0)
        
        plate_primer_long['well_sort_key'] = plate_primer_long['well'].apply(natural_sort_key)
        plate_primer_long = plate_primer_long.sort_values(['well_sort_key', 'name']).reset_index(drop=True)
        plate_primer_long = plate_primer_long.drop(columns=['well_sort_key'])
        
        # Add plate number for identification
        plate_primer_long['plate_number'] = plate_num + 1
        
        print(f"\nPlate {plate_num + 1}:")
        print(plate_primer_long)
        
        all_plates.append(plate_primer_long)
    
    return all_plates


#%%
if __name__ == '__main__':

    # os.chdir('/media/hao/Data/Others/2025-08-15_TF_cloning')
    os.chdir('/mnt/mint/Data/Others/2025-08-15_TF_cloning')

    # by enzyme cut
    res = make_plasmid('ZQ360_Lenti_G11-Xba1-shuttle-BamH1-NLS-Gal4DBD-HT3-T2A-G1-9.dna',
                    insert_dir='test',
                    insert_name_regex='ZQ.*_Lenti_G11-Xba1-(.*)-BamH1-NLS',
                    insert_is_clean=False,
                    enzyme=['XbaI', 'BamHI'],
                    primer_overhang_vector_len=5,
                    primer_overhang_insert_len=25)
    
    # by HR
    res = make_plasmid('ZQ360_Lenti_G11-Xba1-shuttle-BamH1-NLS-Gal4DBD-HT3-T2A-G1-9.dna',
                insert_dir='test',
                insert_name_regex='ZQ.*_Lenti_G11-Xba1-(.*)-BamH1-NLS',
                insert_is_clean=False,
                HR_sequence=['GATCTGGCAGCGGTGGAGGCTCTAGA', 'GGATCCAGTGGCGGGAGTGGAGGTG'],
                primer_overhang_vector_len=25,
                primer_overhang_insert_len=25)

    # by input sequence
    f = '2025-08-29_cloning.csv'
    part1 = pd.read_csv(f) # part1.columns
    res = make_plasmid('ZQ360_Lenti_G11-Xba1-shuttle-BamH1-NLS-Gal4DBD-HT3-T2A-G1-9.dna',
            insert_sequences=part1['cds'].to_list(),
            insert_names=part1['symbol'].to_list(),
            insert_is_clean=True,
            HR_sequence=['GATCTGGCAGCGGTGGAGGCTCTAGA', 'GGATCCAGTGGCGGGAGTGGAGGTG'],
            primer_overhang_vector_len=25,
            primer_overhang_insert_len=28)


    # write snapgene
    out_dir = 'ZQ360'
    os.makedirs(out_dir, exist_ok=True)
    # res.keys()
    for idx, _ in tqdm(enumerate(res['constructs'])):
        write_snapgene(construct=res['constructs'][idx], 
                       path=os.path.join(out_dir, 'ZQ360_Lenti_G11-' + res['names'][idx] + '-NLS-Gal4DBD-HT3-T2A-G1-9' + '.dna'),
                       insert_seq=res['insert_seqs'][idx], 
                       insert_name=res['names'][idx],
                       vector_seq=res['vector_seq'], 
                       vector_feats=res['vector_feats']
                       )
    # write primer
    res['primers'].to_csv(f.replace('.csv','_primer.csv'), index=False)
    # convert to plate IDT format
    primer_plates = convert_to_plate_format(res['primers'])
    for i, plate in enumerate(primer_plates):
        plate.to_csv(f.replace('.csv', f'_primer_plate_{i+1}.csv'), index=False)


    # write snapgene
    res2 = make_plasmid('ZQ361_Lenti_G10-Xba1-shuttle-BamH1-NLS-VPR-HT6.dna',
        insert_sequences=part1['cds'].to_list(),
        insert_names=part1['symbol'].to_list(),
        insert_is_clean=True,
        HR_sequence=['GATCTGGCAGCGGTGGAGGCTCTAGA', 'GGATCCAGTGGCGGGAGTGGAGGTG'],
        primer_overhang_vector_len=25,
        primer_overhang_insert_len=28)
    out_dir2 = 'ZQ361'
    os.makedirs(out_dir2, exist_ok=True)
    for idx, _ in tqdm(enumerate(res2['constructs'])):
        write_snapgene(construct=res2['constructs'][idx], 
            path=os.path.join(out_dir2, 'ZQ361_Lenti_G10-' + res['names'][idx] + '-NLS-VPR-HT6' + '.dna'),
            insert_seq=res2['insert_seqs'][idx], 
            insert_name=res2['names'][idx],
            vector_seq=res2['vector_seq'], 
            vector_feats=res2['vector_feats']
            )





# %%
