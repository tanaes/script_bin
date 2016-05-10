#!/usr/bin/env python
"""
Reads a pair of fastqs, trims each to indicated length, and returns the 
concatenation of the read pairs. e.g.:

in1: ACGCTCATCTAGGGCATACGCGT
     =================------

in2: GCTCTATCGCGATACGCGTATGC
     =================------

out: ACGCTCATCTAGGGCATGCTCTATCGCGATACGC
     1111111111111111122222222222222222

optionally can choose which orientation, e.g. ff, rf, fr, or rr, but NOT 
reverse complement. 

usage:
filter_fastq_by_fastq.py -r1 [fastq f read] -r2 [fastq r read] -t1 150 -t2 150 > out.fastq
"""
from __future__ import print_function
import sys
import gzip
import argparse
from itertools import izip

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-r1', '--read1_fastq', 
    type=str,
    help='path to fastq for read one')

parser.add_argument('-r2', '--read2_fastq', 
    type=str,
    help='path to fastq for read two')

parser.add_argument('-t1', '--read1_trim', 
    type=int,
    help="number of bases from 5' end for read 1 (default all bases)")

parser.add_argument('-t2', '--read2_trim', 
    type=int,
    help="number of bases from 5' end for read 2 (default all bases)")

parser.add_argument('-o', '--orientation', 
    type=str, default='ff',
    help='path to fastq for read two (default: %(default)s')


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


def trim_cat_fastqs(r1_f, r2_f, r1_t, r2_t, orientation='ff'):
    f1 = readfq(r1_f)
    f2 = readfq(r2_f)

    i = 0
    for ((n1,s1,q1),(n2,s2,q2)) in izip(readfq(r1_f),readfq(r2_f)):
        i += 1
        if n1[-2:] in ['/1','/2']:
            n1 = n1[:-2]
        if n2[-2:] in ['/1','/2']:
            n2 = n2[:-2]
        if n1 != n2:
            raise IOError("Sequences don't match! name1: {0} name2: {1}".format(n1,n2))

        if len(s1) < r1_t:
            raise IndexError("Read 1 is too short to trim: "
                              "\n{0}\n{1}\n{2}".format(n1,s1,q1))
        if len(s2) < r2_t:
            raise IndexError("Read 2 is too short to trim: "
                              "\n{0}\n{1}\n{2}".format(n2,s2,q2))

        s1_t = s1[:r1_t]
        q1_t = q1[:r1_t]
        if orientation[0] == 'r':
            s1_t = s1_t[::-1]
            q1_t = q1_t[::-1]

        s2_t = s2[:r2_t]
        q2_t = q2[:r2_t]
        if orientation[1] == 'r':
            s2_t = s2_t[::-1]
            q2_t = q2_t[::-1]

        print('@{0}\n{1}\n+\n{2}'.format(n1, s1_t + s2_t, q1_t + q2_t))

    return(i)


def main():
    args = parser.parse_args()

    read1_fastq = args.read1_fastq
    read2_fastq = args.read2_fastq
    read1_trim = args.read1_trim
    read2_trim = args.read2_trim
    orientation = args.orientation

    try:
        r1_f = gzip.open(read1_fastq)
        r2_f = gzip.open(read2_fastq)
    except IOError:
        try:
            r1_f = open(read1_fastq)
            r2_f = open(read2_fastq)
        except:
            raise IOError('Could not open files')

    i = trim_cat_fastqs(r1_f, r2_f, read1_trim, read2_trim, orientation)

    print('Joined {0} records'.format(i, file=sys.stderr),file=sys.stderr)


if __name__ == "__main__":
    main()
