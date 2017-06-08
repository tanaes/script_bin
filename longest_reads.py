#!/usr/bin/env python
"""
Outputs n longest fastq records. Can take either gzipped or plain fastqs.

usage:
filter_fastq_by_fastq.py -f [fastq] -n 2 > out.fastq
"""
from __future__ import print_function
import sys
import gzip
import argparse
import heapq

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-f', '--fastq', 
    type=str,
    help='path to fastq for read one')

parser.add_argument('-n', '--num_outputs', 
    type=int, default=1,
    help='number of records to output (default: %(default)s')


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


def n_longest(f, n):
    longest = []

    for (n_i,s_i,q_i) in readfq(f):
        if len(longest) < n or len(s_i) > len(longest[0][1]):
            if len(longest) == n: heapq.heappop(longest)
            heapq.heappush(longest, (n_i,s_i,q_i))

    return(longest)


def main():
    args = parser.parse_args()

    fastq = args.fastq
    n = args.num_outputs

    try:
        f = gzip.open(fastq)
        longest = n_longest(f, n)
        print('blah')
    except IOError:
        print('foo')
        try:
            f = open(fastq)
            longest = n_longest(f, n)
        except IOError as e:
            print('Could not open input files:\n%s' % e,file=sys.stderr)
        except Exception as e:
             print(e,file=sys.stderr)
    except Exception as e:
         print(e,file=sys.stderr)

    for record in longest:
        print('@{0}\n{1}\n+\n{2}\n'.format(record[0],record[1],record[2]))

if __name__ == "__main__":
    main()
