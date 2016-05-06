#!/usr/bin/env python
"""
Formats a Humann2 custom database annotation file using the CAZyDB fasta as
downloaded from dbCAN, along with the Ete2 NCBI taxonomy parser and NCBI's
protein accession to taxid file
"""

from ete2 import NCBITaxa
from itertools import groupby
import gzip
import argparse

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-c', '--cazy_fp', 
    type=str,
    help='path to cazy fasta file from dbCAN')

parser.add_argument('-p', '--p2taxid_fp', 
    type=str,
    help='path to prot.accession2taxid.gz from NCBI')

parser.add_argument('-l', '--len_fp', 
    type=str,
    help='path to all.hmm.ps.len from dbCAN')

parser.add_argument('-o', '--output_fp', 
    type=str,
    help='path to output file for writing humann2 custom db annotations')


def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence

    from https://www.biostars.org/p/710/#1412 user brentp
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def read_NCBI_prot_to_taxid_gz(ncbi_fp):
    prot_to_taxid = {}

    with gzip.open(ncbi_fp) as f:
        headers = f.readline()
        for line in f:
            l = line.strip().split('\t')
            prot_to_taxid[l[1]] = l[2]

    return(prot_to_taxid)

def read_NCBI_prot_to_taxid(ncbi_fp):
    prot_to_taxid = {}

    with open(ncbi_fp) as f:
        headers = f.readline()
        for line in f:
            l = line.strip().split('\t')
            prot_to_taxid[l[1]] = l[2]

    return(prot_to_taxid)

def get_taxon_path(taxid, ncbi, ranks=['superkingdom','kingdom','phylum','class','order','family','genus','species'], rank_headers=['sk__','k__','p__','c__','o__','f__','g__','s__']):

    if rank_headers is None:
        rank_headers = ['' for x in ranks]

    tax_path = [x for x in rank_headers]

    lineage = ncbi.get_lineage(taxid)

    lineage_ranks = ncbi.get_rank(lineage)
    lineage_names = ncbi.get_taxid_translator(lineage)

    for level in range(len(lineage)):
        level_id = lineage[level]
        if lineage_ranks[level_id] in ranks:
            tax_path[ranks.index(lineage_ranks[level_id])] += lineage_names[level_id].replace(' ','_')

    return(tax_path)

def read_hmm_len_fp(len_fp):
    len_dict = {}

    with open(len_fp) as f:
        for line in f:
            len_dict[line.strip().split('.')[0]] = int(line.strip().split('\t')[1])*3

    return(len_dict)

def main():
    args = parser.parse_args()

    cazy_fp = args.cazy_fp
    p2taxid_fp = args.p2taxid_fp
    output_fp = args.output_fp
    len_fp = args.len_fp

    ncbi = NCBITaxa()

    # read in ncbi prot to taxid, keep in memory as a dict

    prot_to_taxid = read_NCBI_prot_to_taxid_gz(p2taxid_fp)

    # open 

    hmm_len = read_hmm_len_fp(len_fp)

    # read in cazy

    cazy_f = fasta_iter(cazy_fp)

    # for each cazy,
    with open(output_fp, 'w') as f:
        for header,seq in cazy_f:

            acc = header.split('|')[0]
            fam = header.split('|')[1]

            try:
                gene_length = hmm_len[fam]
            except KeyError:
                gene_length = 1000

            try:
                taxid = prot_to_taxid[acc]
                taxonomy = '.'.join(get_taxon_path(taxid, ncbi))
            except KeyError:
                taxonomy = 'unclassified'

            outline = '{0}\t{1}\t{2}\t{3}\n'.format(header,fam,gene_length,taxonomy)

            f.write(outline)


if __name__ == "__main__":
    main()