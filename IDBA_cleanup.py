#!/usr/bin/env python
"""Script to clean out extra files from IDBA assemblies"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

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

    f_ps = os.listdir(target_dir)

    print(f_ps)

    for f in f_ps:
        if os.path.isfile(f) and not f in files_to_keep:
            print("Deleting %s" % f)
        else:
            print("Keeping %s")



if __name__ == '__main__':
    main()