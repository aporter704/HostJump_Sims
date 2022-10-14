#!/usr/bin/env python
# coding: utf-8

import dendropy as dp
import numpy as np
import re, sys, os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--t', help = 'the nexus tree from MASTER with annotations', required = True)
parser.add_argument('-lh', '--lh', help = 'the location for the human deme', required = True)
parser.add_argument('-pnh', '--pnh', help = 'the probability of sampling non human cases after first human casea', required = False)

args = parser.parse_args()

tree_file_name = args.t
location_human = args.lh
psamp_after_first_human = 0.8

tree = dp.Tree.get_from_path(tree_file_name, 'nexus')

for node in tree.postorder_node_iter():
    if node.is_leaf():
        location =  node.annotations.get_value('location')
        node.taxon.label = node.taxon.label + '_type' + location


# - Get actual number of migration events and their ages
migration_events = []
node_ages = []
for node in tree.postorder_node_iter():
    reaction = node.annotations.get_value('reaction')
    age = node.distance_from_root()
    node_ages.append(age)
    if reaction == 'Migration':
        location = node.annotations.get_value('location')
        migration_events.append([reaction, location, node.child_nodes()[0].annotations.get_value('location'), age])

migration_events = pd.DataFrame(migration_events)
migration_events.columns = ['reaction', 'ancestor_type', 'child_type', 'time']

num_migration_into_human = (migration_events.child_type == location_human).sum()
age_first_migration = migration_events.loc[migration_events.ancestor_type == '0', 'time'].min()
age_of_first_branching = max([i.distance_from_tip() for i in tree.postorder_internal_node_iter() if len(i.child_nodes()) == 2])

# - Export newick tree with metadata in name and no internal nodes
file_name_with_ground_truth = re.sub('[.].+', '', tree_file_name) + '_nMig' + str(num_migration_into_human) + '_migAge' + str(round(age_first_migration, 2)) + '_splitAge' + str(round(age_of_first_branching, 2))

tree.write_to_path(file_name_with_ground_truth + '.newick.tree', 'newick')
print('Newick tree exported to:'+file_name_with_ground_truth + '.newick.tree', 'newick')
print('Number of migration events into humans: ' +  str(num_migration_into_human))

# # Workflow testing complete up to here
# 1. ~~export tree with the data above.~~ 
# 2. prune tips with removed reaction and location 0. These are animals 
# 2.1 Prune tips with removed reaction and location 1. Thesea are humans
# 2.3 prune nodes that are not leaves and which have a single descendants (migrations)
# 2.4 Calculate number of imports, date of first import, and tmrca of sampled tree
# 2.5 Export tree in nexus. This is the 'sampled' tree
# 
# 3 take initial tree and remove only location 1 from the removed reaction
# 3.1 prune all tips with location 0 prior to the first location 1
# 2.5 Export tree in nexus. This is the 'opportunistically sampled' tree
# 
