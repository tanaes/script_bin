#!/usr/bin/env python

import gzip
import sys
import re
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt

def barcode_heatmap(combos):
    all_i5 = set(combos['i5_seq'])
    all_i7 = set(combos['i7_seq'])

    unk_counts = pd.DataFrame(np.empty((len(all_i5), len(all_i7))), index = all_i5, columns = all_i7)

    unk_counts[:] = np.nan

    for i, row in combos.iterrows():
        unk_counts.loc[row['i5_seq'],row['i7_seq']] = row['NumberReads']

    fig, ax = plt.subplots(1, 1, sharex=False, figsize=(60,60))

    sns.heatmap(np.log10(unk_counts),
                ax = ax,
                cmap='jet',
                cbar = False)

    return(fig)


def main():
    fp = sys.argv[1]

    if len(sys.argv) > 2:
        until = sys.argv[2]
    else:
        until = 40000000

    d = {}

    counter = 0
    regex = r'\:([AGCTN]+)\+([AGCTN]+)$'
    with gzip.open(fp, 'rt') as f:
        for line in f:
            counter += 1
            m = re.search(regex, line)

            # print(line)
            if m:
                pair = (m.group(1), m.group(2))
                if pair not in d:
                    d[pair] = 1
                else:
                    d[pair] += 1

            if counter % 1000000 == 0:
                print("parsed %s lines" % counter)

            if counter >= until:
                break

    i5_seq, i7_seq, counts = [], [], []

    for (i7, i5), count in d.items():
        i7_seq.append(i7)
        i5_seq.append(i5)
        counts.append(count)

    combos = pd.DataFrame({'i5_seq': i5_seq,
                           'i7_seq': i7_seq,
                           'NumberReads': counts})

    print('Top 10 i7 sequences:\n\n%s' % combos.groupby('i7_seq').sum().sort_values('NumberReads',
                                                                                  ascending=False).head(10))
    print('Top 10 i5 sequences:\n\n%s' % combos.groupby('i5_seq').sum().sort_values('NumberReads',
                                                                                  ascending=False).head(10))
    print('Top 10 combos:\n\n%s' % combos.sort_values('NumberReads',
                                                    ascending=False).head(10))

    fig = barcode_heatmap(combos)

    fig.savefig('fig.pdf')


if __name__ == "__main__":
    main()