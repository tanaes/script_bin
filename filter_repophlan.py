#!/usr/bin/env python
"""
Program to filter repophlan results by the quality scores.

Example usage:

furcate_tree.py -r repophlan_microbes_wscores.txt
"""

import os
import sys
import argparse
import unittest
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-r', '--repophlan_fp', 
    type=str,
    help='path to repophlan results')

parser.add_argument('-s', '--score_threshold', 
    type=float, default=0.8,
    help='score threshold to include in good file (default: %(default)s)')

parser.add_argument('-o', '--output_fp', 
    type=str, default='./repophlan_filepaths',
    help='path for output (will output [output_fp].good and [output_fp].bad; default: %(default)s)')

parser.add_argument('-c', '--score_cols', 
    type=str, default='score_faa,score_fna,score_rrna,score_trna',
    help='comma delimited list of columns to use for score calculation (default: %(default)s)')

parser.add_argument('-f', '--fp_cols', 
    type=str, default='faa_lname,ffn_lname,fna_lname,frn_lname,score_avg',
    help='comma delimited list of columns to include in output (default: %(default)s)')

parser.add_argument('-t', '--test', 
    action='store_true',
    help='run unittest')


class TestScore(unittest.TestCase):
    def test_calc_avg_score(self):
        input_dict = {'score_faa': {'G000014725': 0.10000000000000001,
          'G000254175': 0.10000000000000001,
          'G000775715': 0.10000000000000001,
          'G000881615': 0.10000000000000001,
          'G000955195': 0.10000000000000001,
          'G001076295': 0.10000000000000001,
          'G001380675': 0.0},
         'score_fna': {'G000014725': 1.0,
          'G000254175': 0.48599999999999999,
          'G000775715': 0.90200000000000002,
          'G000881615': 1.0,
          'G000955195': 1.0,
          'G001076295': 0.78500000000000003,
          'G001380675': 0.0},
         'score_rrna': {'G000014725': 1.0,
          'G000254175': 1.0,
          'G000775715': 1.0,
          'G000881615': 0.0,
          'G000955195': 0.0,
          'G001076295': 0.90000000000000002,
          'G001380675': 0.0},
         'score_trna': {'G000014725': 1.0,
          'G000254175': 0.80000000000000004,
          'G000775715': 0.90000000000000002,
          'G000881615': 0.0,
          'G000955195': 0.0,
          'G001076295': 0.80000000000000004,
          'G001380675': 0.0}}

        exp_scores = {'G000881615': 0.568649417027694, 
            'G000014725': 1, 
            'G000775715': 0.94937522318669, 
            'G000254175': 0.807963889308086, 
            'G001380675': 0, 
            'G000955195': 0.568649417027694, 
            'G001076295': 0.872837563705059}
        
        input_df = pd.DataFrame(input_dict)
        
        obs_scores = calc_avg_score(input_df)
        
        for item in input_df.index:
            self.assertAlmostEqual(obs_scores[item], exp_scores[item])


def calc_avg_score(genome_df, cols = ['score_faa','score_fna','score_rrna','score_trna']):
    sub_df = genome_df[cols].copy()

    sub_df_norm = ((sub_df - sub_df.mean()) / sub_df.std()).mean(axis=1)

    score_avg = (sub_df_norm - sub_df_norm.min()) / (sub_df_norm.max() - sub_df_norm.min())
    
    return(score_avg)


def run_unittests():
    TestScore('test_calc_avg_score').test_calc_avg_score()


def main():

    args = parser.parse_args()

    input_file = args.repophlan_fp
    output_fp = args.output_fp
    score_cols = args.score_cols.strip().split(',')
    fp_cols = args.fp_cols.strip().split(',')
    test = args.test
    score_threshold = args.score_threshold

    if test:
        run_unittests()
        return(0)

    genome_df = pd.read_csv(input_file, sep='\t', header=0, index_col=0)

    genome_df['score_avg'] = calc_avg_score(genome_df, cols=score_cols)

    good_fp = output_fp + '.good'
    bad_fp = output_fp + '.bad'
    genome_df.loc[genome_df['score_avg'] >= score_threshold, fp_cols].to_csv(good_fp, sep = '\t')
    genome_df.loc[genome_df['score_avg'] < score_threshold, fp_cols].to_csv(bad_fp, sep = '\t')


if __name__ == "__main__":
    main()



