#!/usr/bin/env python
"""
Add psuedocount to biom table. Script will add a pseudocount equal to 
[pseudocount scalar] * minimum nonzero entry in biom table. 

[pseudocount scalar] is optional, defaul = 0.5

Use:

add_biom_pseudocount_to_biom.py <input biom fp> <output biom fp> [psuedocount scalar]


"""

import sys
from biom import load_table, Table
import numpy as np

def add_biom_pseudocount(biom_table, pseudo_scalar = 0.5):

    m_dat = biom_table.matrix_data.todense()

    pseudocount = np.min(m_dat[np.nonzero(m_dat)]) * pseudo_scalar

    m_dat[np.where(m_dat == 0)] = pseudocount

    new_table = Table(m_dat, 
                      biom_table.ids(axis='observation'),
                      biom_table.ids(axis='sample'),
                      biom_table.metadata(axis='observation'),
                      biom_table.metadata(axis='sample'),
                      biom_table.table_id)

    return(new_table)


biom_table = load_table(sys.argv[1])

if len(sys.argv) > 3:
    pseudo_scalar = sys.argv[3]
else:
    pseudo_scalar = 0.5

new_table = add_biom_pseudocount(biom_table, pseudo_scalar)

with open(sys.argv[2], 'w') as f:
    f.write(new_table.to_json('Psuedocounts added by add_biom_pseudocount_to_biom.py'))



