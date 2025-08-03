#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 18:22:45 2025

Make plasmid from vector and insert fasta files.

Should run in "ngs" conda environment

@author: hao
"""

import os
from glob import glob
from Bio import SeqIO, Seq
from Bio.Restriction import RestrictionBatch as Batch
import autosnapgene as snap
import re
from tqdm import tqdm



def make_plasmid(work_dir: str,
                 vector_file: str, 
                 insert_dir: str, 
                 replace_left: dict = {'.*[pP][lL]enti_': 'ZQ361-ZQ329 pLenti_GFP10-'},
                 replace_right: dict = {'[-_]EGFP-NLS.*': '-NLS-VPR-HT6.dna'},
                 dest_dirname: str = "new_constructs", 
                 enzyme_list: list[str, str] = ['XbaI', 'BamHI']):
    """
    Make plasmid from vector and insert fasta files.
    Args:
        work_dir: str, path to the work directory
        vector_file: str, path to the vector file
        insert_dir: str, path to the insert directory
        replace_left: dict, left part of the insert name, key is insert left part, value is final left part
        replace_right: dict, right part of the insert name, key is insert right part, value is final right part
        dest_dirname: str, name of the destination directory
        enzyme_list: list[str, str], list of enzymes, only two enzymes are supported
    """

    dest_dirpath = os.path.join(work_dir, dest_dirname)
    os.makedirs(dest_dirpath, exist_ok=True)
    rb = Batch(enzyme_list)

    # vector
    vector = snap.parse(os.path.join(work_dir, vector_file))
    vector_feats = vector.features
    vector_seq = Seq.Seq(vector.sequence)

    # vector cuts
    s1 = rb.get(enzyme_list[0]).search(vector_seq)[0] - 1
    s2 = rb.get(enzyme_list[1]).search(vector_seq)[0] - 1
    vector_cuts = [vector_seq[:s1], # left part
                vector_seq[s2:] # right part
                ]   

    # insert
    insert_files = glob(os.path.join(work_dir, insert_dir, '*.fa'))

    for f in tqdm(insert_files):
        # print(f)
        # f = insert_files[0]
        insert= SeqIO.parse(f, 'fasta')
        for record in insert:
            insert_seq = record.seq

        insert_cuts_1 = rb.get(enzyme_list[0]).search(insert_seq) # left part
        insert_cuts_2 = rb.get(enzyme_list[1]).search(insert_seq) # right part

        if len(insert_cuts_1) == 0 or len(insert_cuts_2) == 0:
            print(f'\n{f} has no cuts')
            continue

        insert_cuts_1 = insert_cuts_1[0] - 1 # left part most left
        insert_cuts_2 = insert_cuts_2[-1] - 1 # right part most right

        insert_seq_cut = insert_seq[insert_cuts_1:insert_cuts_2]

        # final construct sequence
        final_construct = vector_cuts[0] + insert_seq_cut + vector_cuts[1]

        # save name
        insert_basename = os.path.basename(f).replace('.fa', '')
        insert_basename = re.sub(f'{list(replace_left.keys())[0]}', 
                                 '', insert_basename)
        insert_basename = re.sub(f'{list(replace_right.keys())[0]}', 
                                 '', insert_basename)
        final_path = os.path.join(dest_dirpath,
                                f'{list(replace_left.values())[0]}' + 
                                insert_basename + 
                                f'{list(replace_right.values())[0]}')

        # save to snapgene
        tac = snap.SnapGene()
        tac.sequence = str(final_construct)
        tac.set_topology('circular')
        # add vector features
        for feat in vector_feats:
            # feat = vector_feats[2]
            try:
                update_seq = str(vector_seq[feat.segments[0].range[0]-1:feat.segments[0].range[1]])
                if final_construct.find(update_seq):
                    update_feat = snap.Feature.from_segment(
                        name=feat.name, 
                        color=feat.segments[0].color, 
                        # directionality = feat.directionality,
                        type=feat.type)
                    tac.add_feature(update_feat, update_seq)
            except:
                continue

        # add insert features   
        new_feat = snap.Feature.from_segment(
            name=insert_basename, type='CDS', color='#31849b')
        tac.add_feature(new_feat, insert_seq_cut)

        # save
        tac.write(final_path)




if __name__ == '__main__':

    make_plasmid(
        work_dir = '/home/hao/Downloads/construct 21 RTK-KD-VPR (2)',
        vector_file = 'ZQ361-ZQ329.dna',
        insert_dir = 'fasta',
        replace_left = {'.*[pP][lL]enti_': 'ZQ361-ZQ329 pLenti_GFP10-'},
        replace_right = {'[-_]EGFP-NLS.*': '-NLS-VPR-HT6.dna'},
        dest_dirname = "new_constructs", 
        enzyme_list = ['XbaI', 'BamHI']
    )


    make_plasmid(
        work_dir = '/home/hao/Downloads/construct 22 IDR-GAL4DBD (2)',
        vector_file = 'ZQ360-ZQ327.dna',
        insert_dir = 'fasta',
        replace_left = {'.*[pP][lL]enti_': 'ZQ360-ZQ327 pLenti_GFP11-'},
        replace_right = {'[-_]EGFP-NLS.*': '-NLS-Gal4DBD-HT3-T2A-GFP1-9.dna'},
        dest_dirname = "new_constructs", 
        enzyme_list = ['XbaI', 'BamHI']
    )


    make_plasmid(
        work_dir = '/home/hao/Downloads/construct 22 IDR-GAL4DBD (2)',
        vector_file = 'ZQ361-ZQ329.plenti_GFP10-Xba1-shuttle-BamH1-NLS-VPR-HT6.dna',
        insert_dir = 'fasta',
        replace_left = {'.*[pP][lL]enti_': 'ZQ361_pLenti_GFP11-'},
        replace_right = {'[-_]EGFP-NLS.*': '-NLS-VPR-HT6.dna'},
        dest_dirname = "new_constructs", 
        enzyme_list = ['XbaI', 'BamHI']
    )
