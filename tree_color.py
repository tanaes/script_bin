#!/usr/bin/env python

import numpy as np
import math
from colorsys import hls_to_rgb
from StringIO import StringIO
from skbio import TreeNode
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle

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


def render_tree(tree,tip_hues,tip_lums,saturation=0.9):
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

    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False

    ete_tree.show(tree_style = ts)

def layout(node):
    if node.is_leaf():
        # Add node name to laef nodes
        N = AttrFace("name", fsize=14, fgcolor=node.tip_color)
        faces.add_face_to_node(N, node, 0)

ete_tree = Tree()
ete_tree.populate(17,random_branches=True)

sk_tree = TreeNode.read(StringIO(unicode(ete_tree.write())))

tip_hues = get_tip_hues(sk_tree, hue_start=0, hue_end=.8)

#tip_lums = get_tip_lums_by_hue(tip_hues, lum_start=.8, lum_end=.4, hue_per_lum=7)
tip_lums = get_tip_lums_by_clade(sk_tree, depth_thresh=.4)



