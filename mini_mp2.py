#!/usr/bin/env python

import click
import os
import pickle
import bz2
import tempfile
import subprocess
import re
from os.path import join

"""
This will subset the MetaPhlan2 database to only include a subset
of clades you want, to produce a smaller test database.

Requires seqtk and bowtie2 executables in $PATH

Example usage:
==============
mini_mp2.py reduce -b ~/ref_data/metaphlan2/db_v20/mpa_v20_m200 \ 
-p ~/ref_data/metaphlan2/db_v20/mpa_v20_m200.pkl \
-i ~/ref_data/metaphlan2/utils/markers_info.txt.bz2 \
-o test_mini s__Escherichia_coli
"""

@click.group()
def run():
    pass

def _bowtie2_markers(db_path, tmp_dir):
    print('Making temporary fasta output of complete DB')
    markers_fp = join(tmp_dir, 'markers.fasta')

    with open(markers_fp, 'w') as f:
        proc = subprocess.Popen(['bowtie2-inspect', db_path],
                                stdout=f, stderr=subprocess.PIPE)
        output = proc.communicate()

    if proc.returncode != 0:
        raise OSError('Failed with message:\n%s' % output[1])

    return(markers_fp)

def _get_subset_markers(info_fp, clades, tmp_dir):
    print('Getting the subset of markers from desired clades')
    p = re.compile('|'.join(clades))
    mini_markers = []
    # bunzip to temp
    tmp_info_fp = os.path.join(tmp_dir, 'markers_info.txt')
    
    with open(tmp_info_fp, 'w') as f:
        proc = subprocess.Popen(['bunzip2', '-c', info_fp],
                                stdout=f, stderr=subprocess.PIPE)
        output = proc.communicate()

    if proc.returncode != 0:
        raise OSError('Failed with message:\n%s' % output[1])
        
    with open(tmp_info_fp, 'r') as f:
        for line in f:
            if p.search(line):
                mini_markers.append(line.split('\t')[0])

    return(mini_markers)

def _seqtk_filter_markers(markers_fp, mini_markers, tmp_dir):
    print('Filtering marker fasta down to subset markers')
    # write mini markers to file
    minimarker_fp = join(tmp_dir, 'mini_markers.txt')
    with open(minimarker_fp, 'w') as f:
        f.write('\n'.join(mini_markers))

    subset_fasta_fp = join(tmp_dir, 'mini_markers.fasta')

     # filter w seqtk   
    with open(subset_fasta_fp, 'w') as f:
        proc = subprocess.Popen(['seqtk', 'subseq', markers_fp, minimarker_fp],
                                stdout=f, stderr=subprocess.PIPE)
        output = proc.communicate()

    if proc.returncode != 0:
        raise OSError('Failed with message:\n%s' % output[1])

    return(subset_fasta_fp)


def _bowtie2_index(subset_fasta_fp, output_dir, tmp_dir):
    print('Indexing subset marker fasta')
    markers_fp = join(tmp_dir, 'markers.fasta')
    output_base = join(output_dir, os.path.split(output_dir)[1])
    proc = subprocess.Popen(['bowtie2-build', subset_fasta_fp, output_base],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    output = proc.communicate()

    if proc.returncode != 0:
        raise OSError('Failed with message:\n%s' % output[1])

    return(output_base)

def reduce_mp2_pickle(pkl, mini_markers, output_base):
    print('Reducing MetaPhlAn2 db pickle')
    db = pickle.load(bz2.BZ2File(pkl, 'r'))

    # make an initial set of the db with just the explicit markers
    newdb = {}
    newdb['markers'] = {x: db['markers'][x] for x in mini_markers}

    # find the full set of clades associated with these marker taxa
    clades = set()
    for y in newdb['markers']:
        clades.update(set(newdb['markers'][y]['taxon'].split('|')))

    # now iterate through the original db and make sure we have all keys
    # associated with each of these clades
    for x in db['markers']:
        if db['markers'][x]['clade'] in clades:
            newdb['markers'][x] = db['markers'][x]

    # whole taxonomy is small so add full thing for good measure
    newdb['taxonomy'] = db['taxonomy']

    # Save the new mpa_pkl file
    new_pkl_fp = '%s.pkl' % output_base

    out = bz2.BZ2File(new_pkl_fp, 'w')

    pickle.dump(newdb, out, 2)
    out.close()

    return


@run.command()
@click.argument('clades', nargs=-1)
@click.option('--mp2bt2', '-b', required=True, type=click.STRING,
              help='Metaphlan2 bowtie2 database basename')
@click.option('--pkl', '-p', required=True, type=click.Path(exists=True),
              help='Metaphlan2 db pickle corresponding to mp2db')
@click.option('--info_fp', '-i', required=True, type=click.Path(exists=True),
              help='Metaphlan2 db marker info file')
@click.option('--output_dir', '-o', required=True, type=click.STRING,
              help='Directory name for reduced db output')

def reduce(clades, mp2bt2, pkl, info_fp, output_dir):
    import snakemake
    from skbio.io.registry import sniff
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with tempfile.TemporaryDirectory(dir='./') as tmp_dir:
        # make markers output
        markers_fp = _bowtie2_markers(mp2bt2, tmp_dir)

        # get markers subset
        mini_markers = _get_subset_markers(info_fp, clades, tmp_dir)

        # filter fasta to just the subset of markers
        subset_fasta_fp = _seqtk_filter_markers(markers_fp,
                                                mini_markers,
                                                tmp_dir)

        # index with bowtie2
        output_base = _bowtie2_index(subset_fasta_fp, output_dir, tmp_dir)

        # reduce and save new pickle db
        reduce_mp2_pickle(pkl, mini_markers, output_base)

if __name__ == '__main__':
    run()
