#!/usr/bin/env python
"""Script to clean out extra files from IDBA assemblies"""

import os
import sys


def main():

    target_dir = sys.argv[1]

    if len(sys.argv) > 2:
        files_to_keep = sys.argv[2:]
    else:
        files_to_keep = ['contig-100.fa.gz','scaffold.fa.gz','log','begin','end']

    print(target_dir)
    print(files_to_keep)


if __name__ == '__main__':
    main()