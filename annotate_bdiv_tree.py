#!/usr/bin/env python
import os
import argparse
from skbio import TreeNode
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, RectFace, TextFace

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--beta_metrics', 
    type=str,
    help='comma-separated list of beta-diversity measures')

parser.add_argument('-w', '--otu_widths', 
    type=str,
    help='comma-separated list of OTU widths')

parser.add_argument('-i', '--input_dir', 
    type=str,
    default='./',
    help='directory in which bdiv summary folders are found; default=%(default)s')

parser.add_argument('-o', '--output_fp', 
    type=str,
    default='summarized_bdiv_tree.pdf',
    help='Filename for annotated output tree; must end in .pdf, .png, or .svg; default=%(default)s')

parser.add_argument('-t', '--tree_fp', 
    type=str,
    help='Filepath to master tree to annotate; must have nodes named identically to summarized output')


def annotate_tree_with_rugs(tree,results_dict,labels_list,width=10,height=10):

    for node in tree.traverse():
        if not node.is_leaf():
            for row in range(len(labels_list)):
                for col in range(len(labels_list[row])):

                    #print results_dict[node.name][row][col]
                    while node.name.startswith("'") and node.name.endswith("'"):
                        node.name = node.name[1:-1]
                    support =    results_dict[node.name][row][col]
                    if support is not None:
                        color = hls_to_rgb_hex(0,1 - results_dict[node.name][row][col],0)
                    else:
                        color = '#FF0000'

                    node.add_face(RectFace(width, height, fgcolor='#000000', bgcolor=color), column=col, position = "branch-right")

    return


def hls_to_rgb_hex(h,l,s):
    from colorsys import hls_to_rgb
    (r,g,b) = hls_to_rgb(h,l,s)
    hex_rgb = '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))

    return hex_rgb

def newick_coerce(tree_fp):
    f_list = [newick_3, newick_5, add_tip_branches]
    for f in f_list:
        try:
            return f(tree_fp)
        except IOError as (errno, strerror):
            print "I/O error({0}): {1}".format(errno, strerror)
        except:
            continue
    else:
        raise "Error: could not load newick tree file"
        return 0

def newick_3(tree_fp):
    return Tree(tree_fp, format=3)

def newick_5(tree_fp):
    return Tree(tree_fp, format=5)

def add_tip_branches(tree_fp):
    tree = TreeNode.read(tree_fp)

    escape_names = False

    for node in tree.traverse():
        if node.length is None:
            node.length = 0.0
        if not (node.name.startswith("'") and node.name.endswith("'")):
            escape_names = True

    root, ext = os.path.splitext(tree_fp)

    newpath = root + '.relabeled' + ext
    print newpath
    #tree.write(newpath, escape_names=escape_names)
    tree.write(newpath)
    tree = Tree(newpath, format=3)

    return tree

def load_rug_dict(input_dir,beta_metrics,otus):
    
    results_dict = {}
    row = 0
    labels_list = [ [None]*len(otus) for i in range(len(beta_metrics)) ]
    
    for metric in beta_metrics: 
        col = 0 
        for otu in otus:
            support_fp = input_dir + '/' + metric + '_' + str(otu) + '/jackknife_support.txt'
            support_file = open(support_fp,'r')
            
            labels_list[row][col] = metric + '_' + str(otu)
            
            for line in support_file:
                if line.startswith('#'):
                    continue
                
                node,value = line.strip().split('\t')[0:2]
                if node not in results_dict.keys():
                    results_dict[node] = [ [None]*len(otus) for i in range(len(beta_metrics)) ]
                    
                results_dict[node][row][col] = float(value)
                
            col += 1
            
            support_file.close()
        row += 1
    
    return results_dict,labels_list


def main():
    args = parser.parse_args()

    beta_metrics = args.beta_metrics.split(',')
    otu_widths = args.otu_widths.split(',')
    input_dir = args.input_dir
    output_fp = args.output_fp
    tree_fp = args.tree_fp


    nrows = len(beta_metrics)
    ncols = len(otu_widths)


    results_dict, labels_list = load_rug_dict(input_dir, beta_metrics, otu_widths)

    try:
        tree = Tree(tree_fp, format=3)
    except:
        tree = add_tip_branches(tree_fp)

    annotate_tree_with_rugs(tree, results_dict, labels_list)

    ts = TreeStyle()

    for row in range(len(labels_list)):
        for col in range(len(labels_list[row])):
            ts.legend.add_face(TextFace(labels_list[row][col], fsize=20), column=col)

    tree.render(output_fp, tree_style = ts)
    tree.show(tree_style = ts)


if __name__ == "__main__":
    main()

