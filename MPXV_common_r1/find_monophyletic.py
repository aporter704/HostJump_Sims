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

check_lines = open(tree_file_name).readlines()
if len(check_lines) <= 2 and tree_format == 'nexus':
    print('The age of the tree is:')
    print('\'NA\'')
    print('The number of importation events is:')
    print('\'NA\'')
    print('The age of the first importation is:')
    print('\'NA\'')
    print('The first imported sample was collected in:')
    print('\'NA\'')    
    raise SystemExit()

target_type = label
tree = dp.Tree.get_from_path(tree_file_name, tree_format)

for i in tree.leaf_node_iter():
    i.taxon.label = re.sub('^_|^ |_$| $', '', i.taxon.label)

def check_monophyly(node, target_type):
    tip_types = []
    tips = []
    for i in node.leaf_iter():
        tips.append(i)
        tip_types.append(re.split('_| ', i.taxon.label)[1])
    return( [all([j == target_type for j in tip_types]), tips] )

visited_tips = []
importation_nodes = []
for tip in tree.leaf_node_iter():
    if(len(re.findall(target_type, tip.taxon.label)) > 0):
        if(tip in visited_tips):
            continue
        ancestors_of_tip = []
        for ancestor in tip.ancestor_iter():
            is_monophyletic = check_monophyly(ancestor, target_type=target_type)
            ancestors_of_tip.append(ancestor)
            if(not is_monophyletic[0]):
                most_recent_decendants = check_monophyly(ancestors_of_tip[len(ancestors_of_tip)-2], target_type=target_type)[1]
                [visited_tips.append(i) for i in most_recent_decendants]
                importation_nodes.append(ancestor)
                break

tip_ages = [i.distance_from_root() for i in visited_tips]


importation_ages = [i.distance_from_root() for i in importation_nodes]

print('The age of the tree is:')
print(np.max(tip_ages))
print('The number of importation events is:')
print(len(importation_nodes))
print('The age of the first importation is:')
print(np.min(importation_ages))
print('The first imported sample was collected in:')
print(np.min(tip_ages))

