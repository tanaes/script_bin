#!/usr/bin/env python

import sys

from skbio import TreeNode
from biom import parse_table

tree_fp = sys.argv[1]
biom_fp = sys.argv[2]

tree = TreeNode.read(tree_fp)
biom = parse_table(open(biom_fp,'r'))

new_names = {}

for otu,sample in biom.nonzero():
    if otu not in new_names:
        new_names[otu] = '{0}_{1}'.format(otu,sample)
    else:
        new_names[otu] = '{0}_{1}'.format(new_names[otu],sample)

for tip in tree.tips():
    otu = tip.name
    tip.name = new_names[otu]

for node in tree.non_tips():
    node.name = None

print tree


