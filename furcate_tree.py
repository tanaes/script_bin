#!/usr/bin/env python
"""
Program to split the tips of a tree with zero-length branches according to the
values in QIIME mapping file. 

Example usage:

furcate_tree.py -m map.txt -t tree.tre -c 'TreeName,SampleID' -o furcated_tree.tre
"""

import argparse
from cogent import LoadTree
from qiime.parse import parse_mapping_file_to_dict

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--map_fp', 
    type=str,
    help='path to mapping file')

parser.add_argument('-t', '--tree_fp', 
    type=str,
    help='path to tree')

parser.add_argument('-c', '--categories', 
    type=str,
    help='comma delimited list of fields over which to furcate. First field '
         'must match the tip names on the tree. (e.g. TreeName,SampleID'
         'If only one value, will split by SampleID')

parser.add_argument('-o', '--output_fp', 
    type=str,
    help='path to output directory for writing tree images')


def furcate_tree(tree, map_dict, fields):
    if len(fields) == 1:
        fields.append('SampleID')

    for sample in map_dict.keys():
        map_dict[sample]['SampleID'] = sample
        for rank in range(1,len(fields)):
            print "rank: " + str(rank)
            parent_field = fields[rank - 1]
            this_field = fields[rank]
            
            if map_dict[sample][this_field] not in tree.getNodesDict():
                print "inserting field: " + map_dict[sample][this_field] + " in parent: " +  map_dict[sample][parent_field]
                try:
                    tree.getNodesDict()[map_dict[sample][parent_field]].insert(0, map_dict[sample][this_field])
                except KeyError as err:
                    print("Could not find field {0}, value {1} in tree: {2}".format(parent_field,map_dict[sample][parent_field],err))
        #if sample not in tree.getNodesDict():
        #    tree.getNodesDict()[map_dict[sample][fields[-1]]].insert(0,sample)

    return(tree)


def main():
    args = parser.parse_args()

    categories = args.categories
    map_fp = args.map_fp
    tree_fp = args.tree_fp
    output_fp = args.output_fp

    map_dict = parse_mapping_file_to_dict(map_fp)[0]

    fields = categories.split(',')

    tree = LoadTree(tree_fp)

    furcated_tree = furcate_tree(tree, map_dict, fields)

    tree.writeToFile(output_fp)


if __name__ == "__main__":
    main()

