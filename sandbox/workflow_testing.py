#!/usr/bin/env python3

# Requires subsample_trees.py, find_monophyletic.py, treePostProcessing.py, make_lsd_dates.py

import os, sys, re, subprocess

replicates = [i for i in range(5)] # number of replicates to run

template = 'test_template.xml' # Template xml
prefix = 'set_test_template_r' # the prefix we want to use
beast_path = '~/phyloApps/beast2.6.3/bin/beast '
seqgen_path = '~/phyloApps/Seq-Gen-1.3.4/source/seq-gen' # Maybe include clock rate and genome size as arguments here?
iqtree_path = '~/phyloApps/iqtree-2.1.3-Linux/bin/iqtree2' 
lsd_path = '~/phyloApps/lsd-0.3beta/bin/lsd_unix'

for i in replicates[:2]:
    replicate_xml = open(template, 'r').read()
    output_name = prefix+str(i)
    replicate_xml = re.sub('OUTPUT_FILE_NAME', output_name, replicate_xml)
    open(output_name + '.xml', 'w').writelines(replicate_xml)
    os.system(beast_path + output_name+'.xml') # simulate trees
    tree_processing = re.split(' |\n', os.popen('./treePostProcessing.py -t ' + output_name + '.nexus.tree -lh 1 -rc 1e-4 -seqGen ' + seqgen_path).read()) # manipulate trees

    sim_tree_file = tree_processing[-4]
    sim_aln_file = tree_processing[-2]
       
    subsampling_command = './subsample_trees.py -t ' + sim_tree_file + ' -s ' + sim_aln_file + ' -lh '+ '1' + ' -ph ' + '0.8' + ' -pl ' + '0.05' 
    
    subsampled_output = re.split(' |\n', os.popen(subsampling_command).read())

    opp_sampling_aln = subsampled_output[-2]
    opp_sampling_tree = subsampled_output[-3]

    low_samp_aln = subsampled_output[-4]
    low_samp_tree = subsampled_output[-5]

    high_samp_aln = subsampled_output[-6]
    high_samp_tree = subsampled_output[-7]

    estimated_trees = []
    iqtree_command = iqtree_path + ' -s SEQ_ALN_NAME -m GTR+G'
    lsd_command = lsd_path + ' -i SEQ_ALN_NAME -d DATE_ALN_NAME -r a -c'
    for aln in [sim_aln_file, opp_sampling_aln, low_samp_aln, high_samp_aln]:
        os.system(re.sub('SEQ_ALN_NAME', aln, iqtree_command)) # Estimate tree
        os.system('./make_lsd_dates.py -s '+aln)
        os.system(re.sub('DATE_ALN_NAME', re.sub('.fasta', '.date', aln), re.sub('SEQ_ALN_NAME', aln+'.treefile', lsd_command))) # Run lsd
        # capture lsd output and store in variable
        estimated_trees.append(aln+'.treefile.result.nexus')
    
    print(estimated_trees)
    # open data base and store  infor for alles r
        