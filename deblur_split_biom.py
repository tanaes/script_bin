#!/usr/bin/env python
"""
deblur_split_biom.py 
"""

from __future__ import print_function
import sys
import os
import argparse
from biom import load_table

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input_biom', 
    type=str,
    help='path to deblur biom to be split')

parser.add_argument('-f', '--output_fp', 
    type=str,
    help='path stub for output fastqs; will add ".Rx.fastq" (default: biom file basename)')

parser.add_argument('-q1', '--read1_q', 
    type=str, default='I',
    help='quality score(s) to use for read 1 (default: %(default)s)')

parser.add_argument('-q2', '--read2_q', 
    type=str, default='I',
    help='quality score(s) to use for read 2 (default: %(default)s)')

parser.add_argument('-t1', '--read1_trim', 
    type=int,
    help="number of bases from 5' used for read 1 (default 1/2 bases)")

parser.add_argument('-t2', '--read2_trim', 
    type=int,
    help="number of bases from 5' used for read 2 (default 1/2 bases)")

parser.add_argument('-o', '--orientation', 
    type=str, default='ff',
    help='path to fastq for read two (default: %(default)s)')

def uncat_seqs_to_fastq(in_seqs, q1_t, q2_t, r1_t, r2_t, orientation='ff'):
    """
    Generator function for splitting seqs and adding quality scores
    """
    for s in in_seqs:
        if r1_t is None or r2_t is None:
            if len(s) % 2 != 0:
                raise IndexError("Read is not even length, cannot split 1/2")
            r1_t = r2_t = int(len(s) / 2)

        if len(s) < r1_t + r2_t:
            raise IndexError("Read is too short to split to given length: "
                              "\n{0}\n{1}{2}".format(s,'1'*r1_t,'2'*r2_t))

        # set quality score length equal to returned length
        if len(q1_t) != r1_t:
            q1_t = q1_t[0] * r1_t
        if len(q2_t) != r2_t:
            q2_t = q2_t[0] * r2_t

        s1_t = s[:r1_t]
        
        if orientation[0] == 'r':
            s1_t = s1_t[::-1]

        s2_t = s[-r2_t:]
        
        if orientation[1] == 'r':
            s2_t = s2_t[::-1]

        r1 = '@{0}\n{1}\n+\n{2}\n'.format(s, s1_t, q1_t)
        r2 = '@{0}\n{1}\n+\n{2}\n'.format(s, s2_t, q2_t)

        yield r1, r2

def main():
    args = parser.parse_args()

    input_biom = args.input_biom
    output_fp = args.output_fp
    read1_q = args.read1_q
    read2_q = args.read2_q
    read1_trim = args.read1_trim
    read2_trim = args.read2_trim
    orientation = args.orientation

    deblur_biom = load_table(input_biom)

    if output_fp is None:
        output_fp = os.path.splitext(input_biom)[0]

    r1_fp = output_fp + '.R1.fastq'
    r2_fp = output_fp + '.R2.fastq'

    joined_seqs = deblur_biom.ids(axis='observation')
    split_seqs = uncat_seqs_to_fastq(joined_seqs, read1_q, read2_q, read1_trim, read2_trim, orientation=orientation)

    with open(r1_fp, 'w') as r1_f, open(r2_fp, 'w') as r2_f:
        i = 0
        for r1, r2 in split_seqs:
            i += 1
            r1_f.write(r1)
            r2_f.write(r2)

    print('Split {0} records'.format(i, file=sys.stderr),file=sys.stderr)


if __name__ == "__main__":
    main()

