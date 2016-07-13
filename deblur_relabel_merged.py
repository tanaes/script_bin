#!/usr/bin/env python
"""
deblur_relabel_merged.py

This script takes a deblurred biom table with concatenated sequences as IDs
and a merged fastq file that is the result of overlap merging of those
concatenated sequences and relabels the deblurred biom with the new, merged IDs.

Currently this is really slow due to the way the Biom collapses tables -- expect
around 1 minute per 1000 observations.
"""

from __future__ import print_function
import sys
import os
import argparse
from biom import load_table, Table
from biom.parse import biom_open
import unittest

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input_biom_fp', 
    type=str,
    help='path to deblur biom to be split')

parser.add_argument('-o', '--output_biom_fp', 
    type=str,
    help='path to output biom (default: biom file basename + merged.biom)')

parser.add_argument('-f', '--merged_fastq_fp', 
    type=str,
    help='path to merged fastq file (output of deblur_split_biom.py)')

parser.add_argument('-t', '--test', 
    action='store_true',
    help='run unittest')


def readfq(fp): # this is a generator function
    """
    From https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def get_merged_dict(fastq):
    merge_dict = {}

    for name, seq, qual in fastq:
        merge_dict[name] = seq

    return(merge_dict)


def collapse_biom_observations(input_biom, merge_dict):
    output_biom = input_biom.collapse(lambda id_, md: merge_dict[id_],
                                      norm=False, axis='observation')

    return(output_biom)


def main():
    args = parser.parse_args()

    input_biom_fp = args.input_biom_fp
    output_biom_fp = args.output_biom_fp
    merged_fastq_fp = args.merged_fastq_fp
    
    deblur_biom = load_table(input_biom_fp)

    if output_biom_fp is None:
        output_biom_fp = os.path.splitext(input_biom_fp)[0] + '.merged.biom'

    with open(merged_fastq_fp) as fq:
        
        merged_fastq = readfq(fq)

        # read each of the fastqs, make a dict of label:merged read
        merge_dict = get_merged_dict(merged_fastq)

        # filter biom to just the keys of dict
        deblur_biom = deblur_biom.filter(lambda val, id_, md: id_ in merge_dict,
                                         axis='observation')

        output_biom = collapse_biom_observations(deblur_biom, merge_dict)

        with biom_open(output_biom_fp, 'w') as f:
            output_biom.to_hdf5(f, 'deblur_relabel_merged.py')


if __name__ == "__main__":
    main()

