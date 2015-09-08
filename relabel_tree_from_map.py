#!/usr/bin/env python

import argparse
import pandas as pd
from StringIO import StringIO
from skbio import TreeNode

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--mapping_fp', 
    type=str,
    help='path to mapping file')

parser.add_argument('-t', '--tree_fp', 
    type=str,
    help='path to tree file')

parser.add_argument('-f', '--find_field', 
    type=str,
    help='field name in mapping file that matches tips on tree')

parser.add_argument('-r', '--replace_field', 
    type=str,
    help='field name in mapping file to replace tip with')

parser.add_argument('-o', '--output_fp', 
    type=str,
    help='output tree filepath')


def read_mapping_file(map_fp):

    table = pd.read_table(map_fp)

    return(table)


def replace_tips(tree,replace_dict,find_field,replace_field):

    for tip in tree.tips():

        if tip.name in replace_dict:
            print "{0} -> {1}".format(tip.name, replace_dict[tip.name])
            try:

                tip.append(TreeNode.read(StringIO(u"%s;" % replace_dict[tip.name])))
            except:
                tip.name = replace_dict[tip.name]

        else:
            print "%s not in replace_dict" % tip.name

    return tree


def main():
    args = parser.parse_args()
    
    mapping_fp = args.mapping_fp
    tree_fp = args.tree_fp
    find_field = args.find_field
    replace_field = args.replace_field
    output_fp = args.output_fp

    mapping_table = read_mapping_file(mapping_fp)

    replace_dict = dict(zip(mapping_table[find_field],mapping_table[replace_field]))

    tree = TreeNode.read(tree_fp)

    tree = replace_tips(tree,replace_dict,find_field,replace_field)

    tree.write(output_fp)


if __name__ == "__main__":
    main()
