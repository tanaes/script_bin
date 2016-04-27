#!/usr/bin/env python
"""
Reads in pair of fastqs, one of which has been filtered, and returns matching
records for the other. 

usage:
filter_fastq_by_fastq.py -f [filtered fastq] -u [unfiltered fastq] > out.fastq
"""
from __future__ import print_function
import sys
import gzip
import argparse


parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-f', '--filtered_fastq', 
    type=str,
    help='path to fastq that has been filtered')

parser.add_argument('-u', '--unfiltered_fastq', 
    type=str,
    help='path to the fastq that we need to make match')


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


def filter_fastqs(f1_f, f2_f):
    f1 = readfq(f1_f)
    f2 = readfq(f2_f)

    i = 0
    j = 0
    k = 0
    for n1, s1, q1 in f1:
        if n1[-2:] in ['/1','/2']:
            n1 = n1[:-2]
        i += 1
        while True:
            try:
                (n2, s2, q2) = f2.next()
            except StopIteration:
                for n1, s1, q1 in f1:
                    i += 1
                print('Hit end of file!', file=sys.stderr)
                return(i,j,k)
            j += 1
            if n1 == n2:
                print('@{0}\n{1}\n+\n{2}\n'.format(n2, s2, q2))
                k += 1
                break

    return(i,j,k)


def main():
    args = parser.parse_args()

    filtered_fastq = args.filtered_fastq
    unfiltered_fastq = args.unfiltered_fastq

    try:
        f1_f = gzip.open(filtered_fastq)
        f2_f = gzip.open(unfiltered_fastq)

        (i,j,k) = filter_fastqs(f1_f, f2_f)
    except IOError:
        try:
            f1_f = open(filtered_fastq)
            f2_f = open(unfiltered_fastq)

            (i,j,k) = filter_fastqs(f1_f, f2_f)
        except:
            raise 'Could not open file'

    print('Searched through {0} records in {1} for {2} records from {3}, matched {4}'.format(j,unfiltered_fastq,i,filtered_fastq,k), file=sys.stderr)


if __name__ == "__main__":
    main()
