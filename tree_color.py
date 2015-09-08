#!/usr/bin/env python

import numpy as np
import math
import argparse
from colorsys import hls_to_rgb
from StringIO import StringIO
from skbio import TreeNode
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle

parser = argparse.ArgumentParser()

parser.add_argument('-t', '--tree_fp', 
    type=str,
    help='newick-format tree file to color')

parser.add_argument('-o', '--output_fp', 
    type=str,
    default=None,
    help='output image filepath (must have .svg, .pdf, or .png extension')

parser.add_argument('-s', '--supress_display', 
    action='store_true',
    help='supress display of tree in treeview application')

parser.add_argument('--hue_start',
    type=float,
    default=0,
    help='Starting hue, between 0 and 1. default=%(default)s')

parser.add_argument('--hue_end',
    type=float,
    default=.8,
    help='Ending hue, between 0 and 1. default=%(default)s')

parser.add_argument('--saturation',
    type=float,
    default=0.9,
    help='Saturation, between 0 and 1. default=%(default)s')

parser.add_argument('--depth_thresh',
    type=float,
    default=.3,
    help='Clade depth threshold (as a proportion of total tree depth)' \
         ' for iterating lightness values. default=%(default)s')

parser.add_argument('--lum_default',
    type=float,
    default=.7,
    help='Default lighness value for tip labels. default=%(default)s')

parser.add_argument('--lum_max',
    type=float,
    default=.9,
    help='Maximum lightness value for tip labels. default=%(default)s')

parser.add_argument('--lum_min',
    type=float,
    default=.5,
    help='Minimum lightness value for tip labels. default=%(default)s')

parser.add_argument('--lum_sep',
    type=float,
    default=.15,
    help='standard spacing of lightness values for tip labels. default=%(default)s')


def get_tip_hues(tree, hue_start=0, hue_end=.8):
    """Picks hues for tips of tree"""
    # get tip to tip distance matrix
    x = tree.tip_tip_distances()

    # get the postorder tip spacing
    tip_spacing = [x[j,j-1] for j in range(1,x.shape[0])]

    # transform into HSL Hue spacing
    hue_range = hue_end - hue_start

    spacing_hues = np.cumsum(tip_spacing / sum(tip_spacing)) * hue_range
    tip_hues = np.insert(spacing_hues, 0, hue_start)

    return(dict(zip(x.ids, tip_hues)))


def get_tip_lums_by_hue(tip_hues, lum_start=.9, lum_end=.3, hue_per_lum=7):
    """Attempt to pick luminance per tip to maximize separation of similar hues"""

    ### NOTE: this function kind of sucks

    # space out luminance values across the tree, so that trees with few leaves
    # have fewer luminance levels

    num_lum_levels = int(math.ceil(float(len(tip_hues)) / float(hue_per_lum)))

    if num_lum_levels > 1:
        lum_levels = [lum_start - (x * (lum_start - lum_end) / (num_lum_levels - 1)) for x in range(num_lum_levels)]

    else:
        lum_levels = [lum_start]

    # we know how many lum levels there are
    # 

    tip_lums = [0] * len(tip_hues)

    # assign luminosity in alternating pattern, to better separte out neighbor
    # hues.
    l = 0
    while l < len(tip_hues):
        if l < (len(tip_hues) % hue_per_lum) * num_lum_levels:
            tip_lums[l] = lum_levels[l % num_lum_levels]
        else:
            tip_lums[l] = lum_levels[l % (num_lum_levels - 1)]
        l += 1

    return(tip_lums)


def get_tip_lums_by_clade(tree,depth_thresh=.3,lum_default=.7,lum_max=.9,lum_min=.5,lum_sep=.15):
    """Attempt to pick luminance values in clusters"""
    # max distance in tree

    threshold = tree.get_max_distance()[0] * depth_thresh

    tip_names = tree.tip_tip_distances().ids

    tip_group = {}

    # find nodes below a threshold of max distance
    # preorder search
    group_iter = 0
    for node in tree.preorder():
        if node.get_max_distance()[0] < threshold and not node.is_tip():
            for tip in node.tips():
                if tip.name not in tip_group:
                    tip_group[tip.name] = group_iter
            group_iter += 1

    # initialize array of tip lightness
    tip_lums = {tip: lum_default for tip in tip_names}

    # iterate over each tip grouping:
    for group in range(group_iter):
        # get tips for that group
        tips = [k for k, v in tip_group.iteritems() if v == group]
        # if there is space in the defined lightness range to space by default
        # lightness spacing, do that:
        if ((len(tips) - 1) * lum_sep) <= (lum_max - lum_min):
            lum_start = (lum_min + (lum_max - lum_min) / 2) - (((len(tips) - 1) * lum_sep) / 2)
            j = 0
            for tip in tips:
                tip_lums[tip] = lum_start + j * lum_sep
                j += 1
        else:
            lum_step = (lum_max - lum_min) / (len(tips) - 1)
            j = 0
            for tip in tips:
                tip_lums[tip] = lum_min + j * lum_step
                j += 1

    return(tip_lums)


def hls_to_rgb_hex(h,l,s):
    (r,g,b) = hls_to_rgb(h,l,s)
    hex_rgb = '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))

    return hex_rgb


def render_tree(sk_tree, tip_hues, tip_lums, saturation=0.9, output_fp=None, supress_display=False):
    ete_tree = Tree(str(sk_tree))
    tree_tips = [x.name for x in sk_tree.tips()]

    for ete_leaf in ete_tree.iter_leaves():
        if (ete_leaf.name not in tree_tips) and (ete_leaf.name[1:-1] in tree_tips):
            ete_leaf.name = ete_leaf.name[1:-1]
        elif ete_leaf.name not in tree_tips:
            raise 'leaf {0} in ete-parsed tree not found in skbio tree {1}'.format(ete_leaf.name,str(tree))

    for n in ete_tree.traverse():
        if n.is_leaf():
            hex_color = hls_to_rgb_hex(tip_hues[n.name], tip_lums[n.name], saturation)
            n.add_features(tip_color=hex_color)

    style = NodeStyle()
    style["size"] = 0
    for l in ete_tree.traverse():
        l.set_style(style)

    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False

    if output_fp:
        ete_tree.render(output_fp, tree_style = ts)

    if not supress_display:
        ete_tree.show(tree_style = ts)

    return


def layout(node):
    if node.is_leaf():
        # Add node name to leaf nodes
        N = AttrFace("name", fsize=14, fgcolor=node.tip_color)
        faces.add_face_to_node(N, node, 0)


def main():
    args = parser.parse_args()

    tree_fp = args.tree_fp
    output_fp = args.output_fp
    hue_start = args.hue_start
    hue_end = args.hue_end
    saturation = args.saturation
    depth_thresh = args.depth_thresh
    lum_default = args.lum_default
    lum_max = args.lum_max
    lum_min = args.lum_min
    lum_sep = args.lum_sep
    supress_display = args.supress_display

    sk_tree = TreeNode.read(tree_fp)

    tip_hues = get_tip_hues(sk_tree, hue_start=hue_start, hue_end=hue_end)

    #tip_lums = get_tip_lums_by_hue(tip_hues, lum_start=.8, lum_end=.4, hue_per_lum=7)
    tip_lums = get_tip_lums_by_clade(sk_tree,
                                     depth_thresh=depth_thresh,
                                     lum_default=lum_default,
                                     lum_max=lum_max,
                                     lum_min=lum_min,
                                     lum_sep=lum_sep)

    render_tree(sk_tree, tip_hues, tip_lums, saturation=saturation, output_fp=output_fp, supress_display=supress_display)


if __name__ == "__main__":
    main()


