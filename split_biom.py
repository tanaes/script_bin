#!/usr/bin/env python
import argparse
from biom import load_table
from h5py import File
from biom.util import biom_open
from os.path import join

parser = argparse.ArgumentParser()

parser.add_argument("-n", "--n", type=int,
                    help="split biom file into n biom files")
parser.add_argument("-i", "--input_fp",
                    help="biom file to split")
parser.add_argument("-o", "--output_dir",
                    help="dir to put output split biom files")


def main():
    args = parser.parse_args()
    n = args.n
    input_fp = args.input_fp
    outdir = args.output_dir

    biom_table = load_table(input_fp)

    obs_ids = biom_table.ids(axis='observation')

    print "{0} total ids\n".format(len(obs_ids))
    
    chunk_size = int(len(obs_ids)/n)

    last_id = -1

    for chunk in range(1,n):

        begin_id = last_id + 1
        end_id = chunk * chunk_size
        print "chunk: {0} begin: {1} end: {2}\n".format(chunk, begin_id, end_id)

        sub_ids = obs_ids[begin_id : end_id]

        sub_table = biom_table.filter(lambda val, id_, md: id_ in sub_ids, axis='observation', invert=False, inplace=False)
        with biom_open(join(out_dir,'chunk{0}.biom'.format(chunk)), 'w') as out_f:
            sub_table.to_hdf5(out_f, "split_biom.py")

        last_id = end_id

    begin_id = last_id + 1
    chunk += 1

    print "chunk: {0} begin: {1} end: {2}\n".format(chunk, begin_id, len(obs_ids))

    sub_ids = obs_ids[last_id + 1 : ]

    sub_table = biom_table.filter(lambda val, id_, md: id_ in sub_ids, axis='observation', invert=False, inplace=False)
    with biom_open(join(out_dir,'chunk{0}.biom'.format(n)), 'w') as out_f:
        sub_table.to_hdf5(out_f, "split_biom.py")



if __name__ == "__main__":
    main()
