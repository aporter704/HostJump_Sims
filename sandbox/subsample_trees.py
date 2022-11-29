#!/usr/bin/env python
# coding: utf-8

import dendropy as dp
import numpy as np
import re, sys, os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--t', help = 'the nexus tree from MASTER with annotations', required = True)
parser.add_argument('-s', '--s', help = 'the sequence alignment in fasta format', required = True)
parser.add_argument('-lh', '--lh', help = 'location label for human cases', required = True)
parser.add_argument('-ph', '--ph', help = 'high sampling probability', required = True)
parser.add_argument('-pl', '--pl', help = 'low sampling probability', required = True)

args = parser.parse_args()
tree_file_name = args.t
aln_file_name = args.s
location_human = args.lh
non_human_sampling_prop_high = float(args.ph)
non_human_sampling_prop_low = float(args.pl)

# Function to remove samples from alignments
def prune_alignment(aln, labels_to_remove):
    for i in aln:
        if re.sub('\'', '', i.label) in labels_to_remove:
            aln.remove_sequences([i])
    return(aln)

# Read tree and alingment files
print('Reading data')
tree = dp.Tree.get_from_path(tree_file_name, 'nexus')
alignment = dp.DnaCharacterMatrix.get_from_path(aln_file_name, 'fasta')

# Sample tips with high sampling probability and export pruned tree and alignment
print('Sampling non human samples with probability '+str(non_human_sampling_prop_high))
tips_to_remove = []
for i in tree.leaf_node_iter():
    location_temp =  i.annotations.get_value('location')
    if location_temp != location_human:
        to_sample = np.random.binomial(1, non_human_sampling_prop_high)
        if to_sample == 0:
            tips_to_remove.append(i.taxon.label)
            print('Removing '+i.taxon.label)
            
tree_high_sampling = tree.clone(depth = 1)
tree_high_sampling.prune_taxa_with_labels(tips_to_remove)
tree_high_sampling.write_to_path(re.sub('.nexus.tree', '_prunedHighSamp.nexus.tree', tree_file_name), 'nexus')
alignment_high_sampling = alignment.clone(depth = 1)
prune_alignment(alignment_high_sampling, tips_to_remove).write_to_path(re.sub('.fasta', '_prunedHighSamp.fasta', aln_file_name), 'fasta')

# Sample tips with low sampling probility and export pruned tree and alignment
print('Sampling non human samples with probability '+str(non_human_sampling_prop_low))
tips_to_remove = []
for i in tree.leaf_node_iter():
    location_temp =  i.annotations.get_value('location')
    if location_temp != location_human:
        to_sample = np.random.binomial(1, non_human_sampling_prop_low)
        if to_sample == 0:
            tips_to_remove.append(i.taxon.label)
            print('Removing '+i.taxon.label)

tree_low_sampling = tree.clone(depth = 1)
tree_low_sampling.prune_taxa_with_labels(tips_to_remove)
tree_low_sampling.write_to_path(re.sub('.nexus.tree', '_prunedLowSamp.nexus.tree', tree_file_name), 'nexus')
alignment_low_sampling = alignment.clone(depth = 1)
prune_alignment(alignment_low_sampling, tips_to_remove).write_to_path(re.sub('.fasta', '_prunedLowSamp.fasta', aln_file_name), 'fasta')

# Get all tip ages
tip_ages = []
for i in tree.leaf_node_iter():
    tip_ages.append([i.taxon.label, i.annotations.get_value('location'), i.distance_from_root()])

tip_ages = pd.DataFrame(tip_ages)
tip_ages.columns = ['tip_label', 'type', 'date']

# Find age of first human case
first_human_case = tip_ages.date[tip_ages.type == location_human].min()

# Sample tips with high sampling probabilty after first human case
print('Sampling non human samples with probability '+str(non_human_sampling_prop_high)+' after first human case')
tips_to_remove = []
for i in range(tip_ages.shape[0]):
    if tip_ages.date[i] < first_human_case:
        tips_to_remove.append(tip_ages.tip_label[i])
        print('Prunning (opportunistic) '+tip_ages.tip_label[i])
    elif tip_ages.type[i] != location_human and tip_ages.date[i] >= first_human_case:
        if np.random.binomial(1, non_human_sampling_prop_high) == 0:
            tips_to_remove.append(tip_ages.tip_label[i])
            print('Prunning (opportunistic) '+tip_ages.tip_label[i])
            
tree_opportunistic_sampling = tree.clone(depth = 1)
tree_opportunistic_sampling.prune_taxa_with_labels(tips_to_remove)
tree_opportunistic_sampling.write_to_path(re.sub('.nexus.tree', '_prunedOppSamp.nexus.tree', tree_file_name), 'nexus')

aln_opportunistic_sampling = alignment.clone(depth = 1)
prune_alignment(aln_opportunistic_sampling, tips_to_remove).write_to_path(re.sub('.fasta', '_prunedOppSamp.fasta', aln_file_name), 'fasta')
print('Output files')

print(re.sub('.nexus.tree', '_prunedHighSamp.nexus.tree', tree_file_name))
print(re.sub('.fasta', '_prunedHighSamp.fasta', aln_file_name))

print(re.sub('.nexus.tree', '_prunedLowSamp.nexus.tree', tree_file_name))
print(re.sub('.fasta', '_prunedLowSamp.fasta', aln_file_name))

print(re.sub('.nexus.tree', '_prunedOppSamp.nexus.tree', tree_file_name))
print(re.sub('.fasta', '_prunedOppSamp.fasta', aln_file_name))
