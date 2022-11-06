#!/usr/bin/env python
# coding: utf-8

import dendropy as dp
import numpy as np
import re, sys, os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--t', help = 'tree with labels to match. The labels should be the second elementa after an underscore', required = True)
parser.add_argument('-f', '--f', help = 'format of the tree (newick or nexus)', required = True)
parser.add_argument('-l', '--l', help = 'label to match', required = True)

args = parser.parse_args()
tree_file_name = args.t
tree_format = args.f
label = args.l

target_type = label
tree = dp.Tree.get_from_path(tree_file_name, tree_format)

def check_monophyly(node, target_type):
    tip_types = []
    tips = []
    for i in node.leaf_iter():
        tips.append(i)
        tip_types.append(re.split('_', i.taxon.label)[1])
    return( [all([j == target_type for j in tip_types]), tips] )

monophyletic_nodes = []
all_imports = []
visited_tips = []
for tip in tree.leaf_node_iter():
    # if tip in visited_tips:
        continue
    if(re.split('_', tip.taxon.label)[1] == target_type):
        monophyletic_ancestors = []
        for ancestor in tip.ancestor_iter():
            is_monophyletic = check_monophyly(ancestor, target_type) 
            if is_monophyletic[0]:
                monophyletic_ancestors.append(ancestor)
                for child_tip in is_monophyletic[1]:
                    visited_tips.append(child_tip)
            else:
                all_imports.append(ancestor)# Because we want to count singletons too!
                visited_tips.append(tip)
                break
        if len(monophyletic_ancestors) > 0:
            oldest_monophyletic_ancestor = monophyletic_ancestors[-1] 
            if not (oldest_monophyletic_ancestor in monophyletic_nodes):
                monophyletic_nodes.append(oldest_monophyletic_ancestor)
        
print('The number of monophyletic groups including singletons is:')
print(len(all_imports))
print('The number of monophyletic groups excluding singletons is:')
print(len(monophyletic_nodes))