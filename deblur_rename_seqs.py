#!/usr/bin/env python
"""
Reformats a deblurred OTU table with anonymous OTU identifiers (otu1-otuX), 
storing the sequence information as OTU metadata. 

Also optionally outputs a FASTA file of the sequences. 
"""
from __future__ import print_function
import sys
import os
import argparse
import unittest
from biom import Table, load_table
from biom.parse import biom_open
import numpy as np

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input_biom_fp', 
    type=str,
    help='path to deblur biom to be modified')

parser.add_argument('-o', '--output_biom_fp', 
    type=str,
    help='path stub for output fastqs; will add ".renamed.biom" (default: biom file basename)')

parser.add_argument('-f', '--output_fasta_fp', 
    type=str,
    help='path stub for output fastqs (default: none)')

parser.add_argument('-m', '--metadata_name', 
    type=str, default='deblurred_seq',
    help='name to use for sequence metadata category (default: %(default)s)')

parser.add_argument('-n', '--name_stub', 
    type=str, default='deblur',
    help='name to use for sequence headers (default: %(default)s)')

parser.add_argument('-t', '--test', 
    action='store_true',
    help='run unittest')


class TestRename(unittest.TestCase):
    def test_rename_deblur_biom(self):
        input_biom = Table(np.array([[0,0,2],[0,1,3],[3,2,0]]), 
                             ['AAATC','CTTGG','ATCCG'],
                             ['sample1','sample2','sample3'],
                             table_id = 'test')

        exp = Table(np.array([[0,0,2],[0,1,3],[3,2,0]]), 
                             ['deblur0','deblur1','deblur2'],
                             ['sample1','sample2','sample3'],
                             [{'deblurred_seq': 'AAATC'},
                              {'deblurred_seq': 'CTTGG'},
                              {'deblurred_seq': 'ATCCG'}],
                             table_id = 'test renamed')

        obs = rename_deblur_biom(input_biom,
                                         name_stub='deblur',
                                         metadata_name='deblurred_seq')

        self.assertEqual(exp, obs)

    def test_format_seqs_from_deblur_biom(self):
        input_biom_1 = Table(np.array([[0,0,2],[0,1,3],[3,2,0]]), 
                             ['AAATC','CTTGG','ATCCG'],
                             ['sample1','sample2','sample3'],
                             table_id = 'test')

        input_biom_2 = Table(np.array([[0,0,2],[0,1,3],[3,2,0]]), 
                             ['deblur0','deblur1','deblur2'],
                             ['sample1','sample2','sample3'],
                             [{'deblurred_seq': 'AAATC'},
                              {'deblurred_seq': 'CTTGG'},
                              {'deblurred_seq': 'ATCCG'}],
                             table_id = 'test renamed')
        input_biom_3 = Table(np.array([[0,0,2],[0,1,3],[3,2,0]]), 
                             ['deblur0','deblur1','deblur2'],
                             ['sample1','sample2','sample3'],
                             [{'seq': 'AAATC'},
                              {'seq': 'CTTGG'},
                              {'seq': 'ATCCG'}],
                             table_id = 'test renamed')
        obs1 = format_seqs_from_deblur_biom(input_biom_1, name_stub='deblur')
        obs2 = format_seqs_from_deblur_biom(input_biom_2, metadata_name='deblurred_seq')
        exp = '>deblur0\nAAATC\n>deblur1\nCTTGG\n>deblur2\nATCCG'
        self.assertEqual(exp,obs1)
        self.assertEqual(exp,obs2)
        self.assertRaises(TypeError, format_seqs_from_deblur_biom,
                          input_biom_1, metadata_name='deblurred_seq')
        self.assertRaises(KeyError, format_seqs_from_deblur_biom,
                          input_biom_3, metadata_name='deblurred_seq')


def run_unittests():
    TestRename('test_rename_deblur_biom').test_rename_deblur_biom()
    TestRename('test_format_seqs_from_deblur_biom').test_format_seqs_from_deblur_biom()


# def write_anonymous_fasta(seqs, output_fp, name_stub='seq'):
#     with open(output_fp, 'w') as f:
#         i = 0
#         seqnames = [] * len(seqs)
#         for seq in seqs:
#             seqname = name_stub + str(i)
#             seqnames[i] = seqnames
#             f.write('>{0}\n{1}\n'.format(seqname, seq))

#     return(seqnames)


def rename_deblur_biom(biom, name_stub='deblur', metadata_name='deblurred_seq'):
    seqs = biom.ids(axis='observation')

    seqnames = ['{0}{1}'.format(name_stub, x) for x in range(len(seqs))]

    seq_metadata = {seqname: {metadata_name: seq} for seq, seqname in zip(seqs, seqnames)}

    renamed_biom = Table(biom.matrix_data, 
                         seqnames,
                         biom.ids(axis='sample'),
                         biom.metadata(axis='observation'),
                         biom.metadata(axis='sample'),
                         table_id = biom.table_id + ' renamed')

    renamed_biom.add_metadata(seq_metadata, axis='observation')

    return(renamed_biom)


def format_seqs_from_deblur_biom(biom, name_stub='deblur', metadata_name=None):
    if metadata_name is not None:
        try:
            observ_metadata = biom.metadata(axis='observation')
            seqs = [dict(x)[metadata_name] for x in observ_metadata]
            seqnames = biom.ids(axis='observation')
        except KeyError as e:
            #print('Supplied metadata category not present')
            raise e
        except TypeError as e:
            #print('No supplied observation metadata')
            raise e
    else:
        seqs = biom.ids(axis='observation')

        seqnames = ['{0}{1}'.format(name_stub, x) for x in range(len(seqs))]  

    fasta = '\n'.join(['>{0}\n{1}'.format(seqname, seq) for seqname, seq in zip(seqnames, seqs)])
    
    return(fasta)  


def main():
    args = parser.parse_args()

    input_biom_fp = args.input_biom_fp
    output_biom_fp = args.output_biom_fp
    output_fasta_fp = args.output_fasta_fp
    metadata_name = args.metadata_name
    name_stub = args.name_stub
    test = args.test
    
    if test:
        run_unittests()
        return(0)

    deblur_biom = load_table(input_biom_fp)

    if output_biom_fp is None:
        output_biom_fp = os.path.splitext(input_biom_fp)[0] + '.renamed.biom'

    output_biom = rename_deblur_biom(deblur_biom, name_stub=name_stub, metadata_name=metadata_name)

    if output_fasta_fp is not None:
        fasta = format_seqs_from_deblur_biom(output_biom, metadata_name=metadata_name)
        with open(output_fasta_fp, 'w') as f:
            f.write(fasta)
    
    with biom_open(output_biom_fp, 'w') as f:
            output_biom.to_hdf5(f, 'deblur_relabel_merged.py')


if __name__ == "__main__":
    main()
