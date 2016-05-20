#!/usr/bin/env python
"""
Adds full taxonomy strings to Humann2 stratified table.
"""

from ete2 import NCBITaxa
from itertools import groupby
import gzip
import argparse

parser = argparse.ArgumentParser(description=__doc__,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input_fp', 
    type=str,
    help='path to stratifid input Humann2 gene table')

parser.add_argument('-t', '--tdt_out_fp', 
    type=str,
    help='path to output tab-delimited formatted table')

parser.add_argument('-g', '--h2gt_out_fp', 
    type=str,
    help='path to output humann2-formated table')

parser.add_argument('-l', '--rank_headers', 
    type=str, default='sk__,k__,p__,c__,o__,f__,g__,s__',
    help='comma-delimited rank headers to prepend to taxon names'
         'default: %(default)s')

parser.add_argument('-r', '--ranks', 
    type=str, default='superkingdom,kingdom,phylum,class,order,family,genus,species',
    help='comma-delimited ranks to return'
         'default: %(default)s')

def read_humann2_genetable_generator(f):
    header = f.readline().strip().split("\t")
    for l in f:
        line = l.strip().split('\t')
        gene_entry = line[0].split('|')
        gene = gene_entry[0]
        tax = None
        if len(gene_entry) > 1:
            tax = gene_entry[1]
        
        #abund = dict(zip(header[1:],line[1:]))
        yield gene, header[1:], line[1:], tax

def resolve_taxa(higher, lower, ncbi):
    """
    Expects:
    higher: [list] of ncbi taxids
    lower: ncbi taxid
    """

    lineage_lower = ncbi.get_lineage(lower)

    best_higher = higher[0]
    best_score = 0

    for ncbi_id in higher:
        lineage_higher = ncbi.get_lineage(ncbi_id)
        shared = set(lineage_higher) & set(lineage_lower)
        if len(shared) > best_score:
            best_higher = ncbi_id
            best_score = len(shared)
        elif len(shared) == best_score and best_score > 0:
            return(None)

    if best_score:
        return(best_higher)
    else:
        return(None)

def get_best_ncbi_id(family, genus, species, ncbi):
    sp_ncbi = None
    gen_ncbi = None
    fam_ncbi = None

    # get spcies if they have it
    if species is not None:
        sp_ncbi = ncbi.get_name_translator([species])
        if sp_ncbi and len(sp_ncbi[species]) == 1:
            return(sp_ncbi[species][0])
        elif sp_ncbi and len(sp_ncbi[species]) > 1:
            sp_ncbi = sp_ncbi[species]

    # else get genus
    if genus is not None:
        gen_ncbi = ncbi.get_name_translator([genus])
        if gen_ncbi and len(gen_ncbi[genus]) == 1 and not sp_ncbi:
            return(gen_ncbi[genus][0])
        elif gen_ncbi and len(gen_ncbi[genus]) == 1 and len(sp_ncbi) > 1:
            sp_ncbi = resolve_taxa(sp_ncbi, gen_ncbi, ncbi)
            if sp_ncbi:
                return(sp_ncbi)
            else:
                return(gen_ncbi)
        elif gen_ncbi and len(gen_ncbi[genus]) > 1:
            gen_ncbi = gen_ncbi[genus]        

    # else get family
    if family is not None:
        fam_ncbi = ncbi.get_name_translator([family])
        if fam_ncbi and len(fam_ncbi[family]) == 1 and not gen_ncbi:
            return(fam_ncbi[family][0])
        elif fam_ncbi and len(fam_ncbi[family]) == 1 and len(gen_ncbi) > 1:
            gen_ncbi = resolve_taxa(gen_ncbi, fam_ncbi, ncbi)
            if gen_ncbi:
                return(gen_ncbi)
            else:
                return(fam_ncbi)
        elif fam_ncbi and len(fam_ncbi[family]) > 1:
            fam_ncbi = fam_ncbi[family]        
    
    # if can't resolve with nearest higher-level taxon, resolve with bacteria:
    if sp_ncbi and len(sp_ncbi) > 1:
        sp_ncbi = resolve_taxa(sp_ncbi, 2, ncbi)
        if sp_ncbi:
            return(sp_ncbi)

    if gen_ncbi and len(gen_ncbi) > 1:
        gen_ncbi = resolve_taxa(gen_ncbi, 2, ncbi)
        if gen_ncbi:
            return(gen_ncbi)

    if fam_ncbi and len(fam_ncbi) > 1:
        fam_ncbi = resolve_taxa(fam_ncbi, 2, ncbi)
        if fam_ncbi:
            return(fam_ncbi)
    return


def clean_humann2_taxon(tax):
    taxlist = tax.split('.')

    family = None
    if tax == 'unclassified':
        return(None,None,None)
    #print(tax)

    genus = taxlist[0].split('__')[1].split('_')
    #print(genus)
    if len(genus) > 1 and genus[1] == 'noname':
        family = genus[0]
        genus = None
    else:
        genus = genus[0]
  
    species = taxlist[1].split('__')[1].split('_')
    #print(species)
    if len(species) > 2 and species[1] == 'bacterium':
        species = None
    else:
        species = ' '.join(species[:2])
        
    return(family, genus, species)

def get_taxon_path(taxid, ncbi, ranks=['superkingdom','kingdom','phylum','class','order','family','genus','species'], rank_headers=['sk__','k__','p__','c__','o__','f__','g__','s__']):

    if rank_headers is None:
        rank_headers = ['' for x in ranks]

    tax_path = [x for x in rank_headers]

    lineage = ncbi.get_lineage(taxid)

    lineage_ranks = ncbi.get_rank(lineage)
    lineage_names = ncbi.get_taxid_translator(lineage)

    for level in range(len(lineage)):
        level_id = lineage[level]
        if lineage_ranks[level_id] in ranks:
            tax_path[ranks.index(lineage_ranks[level_id])] += lineage_names[level_id].replace(' ','_')

    return(tax_path)


def main():
    args = parser.parse_args()

    input_fp = args.input_fp
    tdt_out_fp = args.tdt_out_fp
    h2gt_out_fp = args.h2gt_out_fp
    rank_headers = args.rank_headers
    ranks = args.ranks

    ncbi = NCBITaxa()

    # input_fp = './R1_trimmed_CAT_rare_genefamilies_cpm_ko.tsv'

    h2gt = read_humann2_genetable_generator(open(input_fp))
    
    rank_headers = rank_headers.split(',')
    ranks = ranks.split(',')
    
    tax_dict = {}

    if tdt_out_fp:
        tdt_out_f = open(tdt_out_fp, 'w')
    if h2gt_out_fp:
        h2gt_out_f = open(h2gt_out_fp, 'w')        

    first_h = True
    first_t = True
    for gene, header, line, tax in h2gt:
        lineage = list(rank_headers)

        if tax and tax not in tax_dict:
            best_id = None
            
            family, genus, species = clean_humann2_taxon(tax)
            best_id = get_best_ncbi_id(family, genus, species, ncbi)

            if best_id is not None:
                lineage = get_taxon_path(best_id, ncbi, ranks=ranks, rank_headers=rank_headers)
            
            tax_dict[tax] = lineage

        if tax:
            lineage = tax_dict[tax]

        if tdt_out_fp:
            if first_t:
                first_t = False
                tdt_out_f.write('Gene Family\t{0}\t{1}\n'.format('\t'.join(header),'\t'.join(ranks)))
            elif tax:
                tdt_out_f.write('{0}\t{1}\t{2}\n'.format(gene,'\t'.join(line),'\t'.join(lineage)))

        if h2gt_out_fp:
            if first_h:
                h2gt_out_f.write('# Gene Family\t{0}\n'.format('\t'.join(header)))
                first_h = False
            if tax:
                if lineage == rank_headers:
                    h2gt_out_f.write('{0}|{1}\t{2}\n'.format(gene,'unknown','\t'.join(line)))
                else:
                    h2gt_out_f.write('{0}|{1}\t{2}\n'.format(gene,'.'.join(lineage),'\t'.join(line)))
            else:
                h2gt_out_f.write('{0}\t{1}\n'.format(gene,'\t'.join(line)))

        

    if tdt_out_fp:
        tdt_out_f.close()
    if h2gt_out_f:
        h2gt_out_f.close()   


if __name__ == "__main__":
    main()