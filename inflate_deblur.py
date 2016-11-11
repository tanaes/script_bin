#!/usr/bin/env python
from __future__ import print_function

import sys
from biom import load_table

deblur_biom = load_table(sys.argv[1])

sample_counts = dict(zip(deblur_biom.ids(axis='sample'),[0] * len(deblur_biom.ids(axis='sample'))))

for seq,sample in deblur_biom.nonzero():
    for i in range(int(deblur_biom.get_value_by_ids(seq,sample))):
        sample_counts[sample] += 1
        print(">{0}_{1}\n{2}".format(sample,sample_counts[sample],seq))
