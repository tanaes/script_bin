#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import gzip


parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-f', '--fastq',
    required=True,
    type=str,
    help='path to fastq')

parser.add_argument('-m', '--min_len', 
    type=int, default=1000,
    help='minimum length to output (default: %(default)s bp')

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

def write_filtered(f):
    for n,s,q in readfq(f):
        if len(s) >= min_len:
            print('@{0}\n{1}\n+\n{2}\n'.format(n,s,q))

if __name__ == '__main__':

    args = parser.parse_args()

    fastq = args.fastq
    min_len = args.min_len


    try:
        f = gzip.open(fastq)
        write_filtered(f)
    except IOError:
        try:
            f = open(fastq)
            write_filtered(f)
        except IOError as e:
            print('Could not open input files:\n%s' % e,file=sys.stderr)
        except Exception as e:
             print(e,file=sys.stderr)
    except Exception as e:
         print(e,file=sys.stderr)
