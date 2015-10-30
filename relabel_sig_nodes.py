#!/usr/bin/env python

import os
import pandas as pd
import collections
from StringIO import StringIO
from biom import load_table
from skbio import TreeNode
import argparse 

# re-parse codiv output file to add cOTU taxonomy info, for phylosift 

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--results_table_fp', 
    type=str,
    help='path to summarized results table to illustrate')

parser.add_argument('-o', '--output_fp', 
    type=str,
    help='path to output results table fp')

parser.add_argument('-s', '--subcluster_dir', 
    type=str,
    help='path to subclustered OTU directory')

parser.add_argument('-p', '--pct', 
    type=float, default=0.5,
    help='percentile for taxonomy resolution (default: %(default)s)')

parser.add_argument('-c', '--collapse_field', 
    type=str,
    default=None,
    help='field name if using collapsed cOTU tables (default: %(default)s)')

parser.add_argument('--top_clade', 
    action='store_true',
    help='only output the deepest node for a series of significant nodes')

def read_cosp_nodes_table(table_fp):

    table = pd.read_table(table_fp)

    return(table)

def resolve_taxonomy(taxonomy_list, pct = 0.5):
    levels = max([len(x) for x in taxonomy_list])

    taxonomy = ['NA'] * levels
    for i in range(levels):
        rank_i_set = [x[i] for x in taxonomy_list]
        
        if collections.Counter(rank_i_set).most_common()[0][1] >= pct * len(rank_i_set):
            taxonomy[i] = collections.Counter(rank_i_set).most_common()[0][0]

    return taxonomy

def main():
    args = parser.parse_args()

    collapse_field = args.collapse_field
    results_table_fp = args.results_table_fp
    subcluster_dir = args.subcluster_dir
    output_fp = args.output_fp
    top_clade = args.top_clade
    pct = args.pct
        
    # read in nodes output file
    results_table = read_cosp_nodes_table(results_table_fp)

    # read in pOTU, load cOTU biom file (for taxonomy)

    # get tips from the s_nodes 

    # get best taxonomy available

    # add convenience check to only print biggest parent node

    # print to output file

    if collapse_field:
        cotu_biom_filename = 'otu_table_%s.biom' % collapse_field
    else:
        cotu_biom_filename = 'otu_table.biom'

    old_tips = []
    old_potu = ''

    i = 0

    for node in results_table.iterrows():

        cotu_tree = TreeNode.read(StringIO(unicode(node[1].s_nodes)), convert_underscores=False)

        tips = [x.name for x in cotu_tree.tips()]

        if set(tips).issubset(old_tips) and top_clade:
            results_table.drop(node[0], inplace=True)
            print 'Dropping row %s' % node[0]
            continue

        if node[1].pOTU is not old_potu:
            cotu_table = load_table(os.path.join(subcluster_dir, str(node[1].pOTU), cotu_biom_filename))

        taxonomy_list = [cotu_table.metadata(x, axis='observation')['taxonomy'] for x in tips]
        
        taxonomy = resolve_taxonomy(taxonomy_list, pct=pct)

        results_table.ix[node[0],'taxonomy'] = '; '.join(taxonomy)

        i += 1

        old_tips = tips

    results_table.to_csv(path_or_buf=output_fp, sep='\t', index=False)

if __name__ == "__main__":
    main()

