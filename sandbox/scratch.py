import os, sys, re, subprocess


subsampling_command = './subsample_trees.py -t set_test_template_r1_nMig4_migAge0.65_splitAge1.98.nexus.tree -s set_test_template_r1_nMig4_migAge0.65_splitAge1.98.fasta -lh 1 -ph 0.8 -pl 0.05'

subsampled_output = re.split('\n', os.popen(subsampling_command).read())

print(subsampled_output[-2])
print(subsampled_output[-3])

print(subsampled_output[-4])
print(subsampled_output[-5])

print(subsampled_output[-6])
print(subsampled_output[-7])







# Get names of trees and alignments