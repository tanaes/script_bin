#!/usr/bin/env python
"""
Reads a fastq split by deblur_fastq_trim_cat.py and splits it back:

in:   ACGCTCATCTAGGGCATGCTCTATCGCGATACGC
      1111111111111111122222222222222222

out1: ACGCTCATCTAGGGCATACGCGT
      =================------

out2: GCTCTATCGCGATACGCGTATGC
      =================------

usage:
deblur_fastq_uncat.py -i [catted fastq] -r1 [fastq f out read] -r2 [fastq r out read] -t1 150 -t2 150
"""
from __future__ import print_function
import sys
import gzip
import argparse
from itertools import izip

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input_fastq', 
    type=str,
    help='path to fastq for input')

parser.add_argument('-r1', '--read1_fastq', 
    type=str,
    help='path to write fastq for read one')

parser.add_argument('-r2', '--read2_fastq', 
    type=str,
    help='path to write fastq for read two')

parser.add_argument('-t1', '--read1_trim', 
    type=int,
    help="number of bases from 5' used for read 1 (default 1/2 bases)")

parser.add_argument('-t2', '--read2_trim', 
    type=int,
    help="number of bases from 5' used for read 2 (default 1/2 bases)")

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


def uncat_fastqs(in_f, r1_f, r2_f, r1_t, r2_t, orientation='ff'):
    i = 0
    f = readfq(in_f)
    for (n,s,q) in readfq(in_f):
        i += 1

        if r1_t is None or r2_t is None:
            if len(s) % 2 != 0:
                raise IndexError("Read is not even length, cannot split 1/2")
            r1_t = r2_t = int(len(s) / 2)

        if len(s) < r1_t + r2_t:
            raise IndexError("Read is too short to split to given length: "
                              "\n{0}\n{1}\n{2}\n{3}{4}".format(n,s,q,'1'*r1_t,'2'*r2_t))

        s1_t = s[:r1_t]
        q1_t = q[:r1_t]
        if orientation[0] == 'r':
            s1_t = s1_t[::-1]
            q1_t = q1_t[::-1]

        s2_t = s[-r2_t:]
        q2_t = q[-r2_t:]
        if orientation[1] == 'r':
            s2_t = s2_t[::-1]
            q2_t = q2_t[::-1]

        r1_f.write('@{0}\n{1}\n+\n{2}\n'.format(n, s1_t, q1_t))
        r2_f.write('@{0}\n{1}\n+\n{2}\n'.format(n, s2_t, q2_t))

    return(i)


def main():
    args = parser.parse_args()

    input_fastq = args.input_fastq
    read1_fastq = args.read1_fastq
    read2_fastq = args.read2_fastq
    read1_trim = args.read1_trim
    read2_trim = args.read2_trim
    orientation = args.orientation

    i = 0

    try:
        in_f = gzip.open(input_fastq)
        r1_f = gzip.open(read1_fastq, 'wb')
        r2_f = gzip.open(read2_fastq, 'wb')
        i = uncat_fastqs(in_f, r1_f, r2_f, read1_trim, read2_trim, orientation)
        in_f.close()
        r1_f.close()
        r2_f.close()
    except IOError:
        try:
            r1_f.close()
            r2_f.close()
            in_f = open(input_fastq)
            r1_f = open(read1_fastq, 'w')
            r2_f = open(read2_fastq, 'w')
            i = uncat_fastqs(in_f, r1_f, r2_f, read1_trim, read2_trim, orientation)
            in_f.close()
            r1_f.close()
            r2_f.close()
        except IOError as e:
            print('Could not open input files:\n%s' % e,file=sys.stderr)
        except Exception as e:
             print(e,file=sys.stderr)
    except Exception as e:
         print(e,file=sys.stderr)

    print('Split {0} records'.format(i, file=sys.stderr),file=sys.stderr)


if __name__ == "__main__":
    main()
