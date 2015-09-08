#!/usr/bin/env python
import argparse
import os
import pandas as pd
from skbio import TreeNode
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace
from biom import parse_table, load_table


parser = argparse.ArgumentParser()

parser.add_argument('-p', '--pies', 
    action='store_true',
    help='annotate cOTU tree with pie charts colored by host')

parser.add_argument('-l', '--labels', 
    action='store_true',
    help='annotate cOTU tree with tip labels colored by host')

parser.add_argument('-c', '--collapse_field', 
    type=str,
    default=None,
    help='field name if using collapsed cOTU tables')

parser.add_argument('-r', '--results_table_fp', 
    type=str,
    help='path to summarized results table to illustrate')

parser.add_argument('-s', '--subcluster_dir', 
    type=str,
    help='path to subclustered OTU directory')

parser.add_argument('-t', '--host_tree_fp', 
    type=str,
    help='path to newick-format host tree for color info')

parser.add_argument('-o', '--output_dir', 
    type=str,
    help='path to output directory for writing tree images')

parser.add_argument('--potu_table_fp', 
    type=str,
    help='path to pOTU table for adding taxonomy information to cOTU trees')

parser.add_argument('-f', '--output_format', 
    choices=['pdf', 'png', 'svg'],
    default='pdf',
    help='format for tree image output. default=%(default)s')

parser.add_argument('--force', 
    action='store_true',
    help='force overwrite of output directory')
# read symbiont trees output

def read_cosp_nodes_table(table_fp):

    table = pd.read_table(table_fp)

    return(table)


# read in mapping file for annotating tips with host data

def read_mapping_file(map_fp):

    table = pd.read_table(map_fp)

    return(table)


# read in host tree + calculate colors
    
def get_host_colors_from_tree(host_tree_fp,
                              hue_start=0,
                              hue_end=.8,
                              depth_thresh=.3,
                              lum_default=.7,
                              lum_max=.9,
                              lum_min=.5,
                              lum_sep=.15,
                              sat=.8):

    from tree_color import get_tip_hues, get_tip_lums_by_clade, hls_to_rgb_hex

    host_tree = TreeNode.read(host_tree_fp)

    host_hues = get_tip_hues(host_tree, hue_start=hue_start, hue_end=hue_end)

    host_lums = get_tip_lums_by_clade(host_tree,
                                      depth_thresh=depth_thresh,
                                      lum_default=lum_default,
                                      lum_max=lum_max,
                                      lum_min=lum_min,
                                      lum_sep=lum_sep)

    host_colors = {x: hls_to_rgb_hex(host_hues[x], host_lums[x], sat) for x in host_hues}

    return(host_colors)

def add_colored_host_label(cotu_tree, cotu_table, host_colors, count = True):
    for tip in cotu_tree.iter_leaves():
        # filter cotu_table to this cotu
        # iterate across nonzero entries and add a host face for that one

        i = 1
        
        for thing in cotu_table.filter(lambda val, id_, md: id_ == tip.name, 
                                       axis='observation',
                                       inplace=False).nonzero():
            host_lab = thing[1]
            
            if count:
                host_lab += "_%s" % cotu_table.get_value_by_ids(thing[0],thing[1])

            tip.add_face(TextFace(host_lab, fgcolor = host_colors[thing[1]]), column=i, position = "branch-right")

            #i += 1

    return


def add_host_pies(cotu_tree, cotu_table, host_colors, count = True, max_pie=100):

    size_max = max(cotu_table.sum(axis='observation'))

    for tip in cotu_tree.iter_leaves():
        # filter cotu_table to this cotu
        # iterate across nonzero entries and add a host face for that one

        pie_colors = [host_colors[host] for host in cotu_table.ids(axis='sample')] 

        thing = cotu_table.filter(lambda val, id_, md: id_ == tip.name, 
                                   axis='observation',
                                   inplace=False)

        size = (((thing.sum(axis='observation')[0]/float(size_max))/3.14) ** (0.5)) * max_pie

        pie_values = [thing.norm(axis='observation').get_value_by_ids(tip.name,host) * 100 for host in thing.ids(axis='sample')] 

        pie = faces.PieChartFace(pie_values,
                              colors=pie_colors,
                              width=size, height=size)

        pie.border.width = None
        pie.opacity = 0.8

        tip.add_face(pie, 1, position="aligned")

    return

def add_taxonomy_label(cotu_tree, potu_table):

    ts = TreeStyle()

    for row in range(len(labels_list)):
        for col in range(len(labels_list[row])):
            ts.legend.add_face(TextFace(labels_list[row][col], fsize=20), column=col)


def print_colored_host_tree(host_tree_fp, host_colors, output_fp):
    tree = Tree(host_tree_fp)

    for host in host_colors:
        if tree.get_leaves_by_name(host) == []:
            tip = "'%s'" % host
        else:
            tip = host

        node = tree.get_leaves_by_name(tip)[0]

        node.add_face(AttrFace("name", fsize=14, fgcolor=host_colors[host]), 0)
    
    # remove stupid blue circles on nodes
    style = NodeStyle()
    style["size"] = 0
    for l in tree.traverse():
        l.set_style(style)

    # get rid of default rendered leaf labels
    ts = TreeStyle()
    ts.show_leaf_name = False
    
    tree.render(output_fp, tree_style = ts)

    return
        

def remove_newick_node_labels(tree_string):
    import re

    return(re.sub(r'\)[\w.]+?\:',r'):', tree_string))


def render_cotu_tree(cotu_tree, output_fp, ts=TreeStyle()):
    ts.show_leaf_name = True
    # ts.scale =  300 # 120 pixels per branch length unit
    # ts.mode = "c"

    nstyle = NodeStyle()
    nstyle["shape"] = "circle"
    nstyle["size"] = 0
    for n in cotu_tree.traverse():
       n.set_style(nstyle)

    print 'rendering tree'
    
    cotu_tree.render(output_fp, tree_style=ts)
# read in otu table for annotating tips with abundance info

# annote symbiont tree by host 

def main():
    args = parser.parse_args()

    collapse_field = args.collapse_field
    results_table_fp = args.results_table_fp
    host_tree_fp = args.host_tree_fp
    potu_table_fp = args.potu_table_fp
    subcluster_dir = args.subcluster_dir
    output_dir = args.output_dir
    pies = args.pies
    labels = args.labels
    output_format = args.output_format
    force = args.force

    if potu_table_fp is not None:
        taxonomy = True
        potu_table = load_table(potu_table_fp)

    results_table = read_cosp_nodes_table(results_table_fp)

    host_colors = get_host_colors_from_tree(host_tree_fp)

    try:
        os.makedirs(output_dir)
    except OSError:
        if not force:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            raise OSError("Output directory already exists. Please choose "
                "a different directory, or force overwrite with -f.")

    host_tree_output_fp = os.path.join(output_dir,'{0}.{1}'.format(
                    os.path.splitext(
                                     os.path.basename(host_tree_fp)
                                    )[0],
                    output_format))

    print_colored_host_tree(host_tree_fp, host_colors, host_tree_output_fp)

    if collapse_field:
        cotu_biom_filename = 'otu_table_%s.biom' % collapse_field
    else:
        cotu_biom_filename = 'otu_table.biom'

    for potu in results_table.iterrows():
        cotu_tree = Tree(remove_newick_node_labels(potu[1].s_nodes))

        cotu_table = load_table(os.path.join(subcluster_dir, str(potu[1].pOTU), cotu_biom_filename))

        print 'Adding hosts to pOTU %s' % potu[1].pOTU

        if pies:
            add_host_pies(cotu_tree, cotu_table, host_colors)

        if labels:
            add_colored_host_label(cotu_tree, cotu_table, host_colors)

        ts = TreeStyle()

        if taxonomy:
            tax = potu_table.metadata(id=str(potu[1].pOTU), axis='observation')['taxonomy']
            tax_string = '\n'.join(tax)
            ts.legend.add_face(TextFace(tax_string), column=0)

        output_fp = os.path.join(output_dir,'{0}.{1}'.format(str(potu[1].pOTU),output_format))

        render_cotu_tree(cotu_tree, output_fp, ts)


if __name__ == "__main__":
    main()




