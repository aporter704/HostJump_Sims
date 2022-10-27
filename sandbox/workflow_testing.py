#!/usr/bin/env python3

import os, sys, re, subprocess

replicates = [i for i in range(5)] # number of replicates to run

template = 'MPXV_speciesjump_common.xml' # Template xml
prefix = 'set1_MPXV_common_r' # the prefix we want to use
beast_path = '~/phyloApps/beast2.6.3/bin/beast '
seqgen_path = '~/phyloApps/Seq-Gen-1.3.4/source/seq-gen'
# Maybe include clock rate and genome size here?

for i in replicates[:4]:
    replicate_xml = open(template, 'r').read()
    output_name = prefix+str(i)
    replicate_xml = re.sub('OUTPUT_FILE_NAME', output_name, replicate_xml)
    open(output_name + '.xml', 'w').writelines(replicate_xml)
    os.system(beast_path + output_name+'.xml')
    tree_processing = re.split(' |\n', os.popen('./treePostProcessing.py -t ' + output_name + '.nexus.tree -lh 1 -rc 1e-4 -seqGen ' + seqgen_path).read())
    sim_tree_file = tree_processing[-3]
    sim_aln_file = tree_processing[-2]
    print(sim_tree_file)
    # Run subsample_trees.py which is ready now
# We need a script to estimate trees
# We need a script to calculate importation events  


#for i in $replicaes; do
#    sed -i.bak 's/OUTPUT_FILE_NAME/$i/g' MPXV_speciesjump_common.xml
 #   done