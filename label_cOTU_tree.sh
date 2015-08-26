#!/bin/bash

FastTree -nt subclustered_otus/${1}/muscle_aligned_seqs/seqs_rep_set_aligned.fasta > subclustered_otus/${1}/${1}_rep_set_aligned.tre

label_tree_by_biom.py subclustered_otus/${1}/${1}_rep_set_aligned.tre subclustered_otus/${1}/otu_table_Species_tree.biom > subclustered_otus/${1}/${1}_rep_set_aligned_labeled.tre


