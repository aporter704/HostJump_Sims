#!/usr/bin/env python3

import os, sys, re, subprocess

replicates = [i for i in range(5)] # number of replicates to run

template = 'test_template.xml' # Template xml
prefix = 'set_test_template_r' # the prefix we want to use
beast_path = '~/phyloApps/beast2.6.3/bin/beast '
seqgen_path = '~/phyloApps/Seq-Gen-1.3.4/source/seq-gen' # Maybe include clock rate and genome size as arguments here?
iqtree_path = ''#  

for i in replicates[:2]:
    replicate_xml = open(template, 'r').read()
    output_name = prefix+str(i)
    replicate_xml = re.sub('OUTPUT_FILE_NAME', output_name, replicate_xml)
    open(output_name + '.xml', 'w').writelines(replicate_xml)
    os.system(beast_path + output_name+'.xml')
    tree_processing = re.split(' |\n', os.popen('./treePostProcessing.py -t ' + output_name + '.nexus.tree -lh 1 -rc 1e-4 -seqGen ' + seqgen_path).read())
    sim_tree_file = tree_processing[-4]
    sim_aln_file = tree_processing[-2]       
    subsampling_command = './subsample_trees.py -t ' + sim_tree_file + ' -s ' + sim_aln_file + ' -lh '+ 'type1' + ' -ph ' + '0.05' + ' -pl ' + '0.05' 
    os.system(subsampling_command)
    ## I need to check the output of the subsampling command!!!
    #subsample_trees.py -t ANNOTATED_NEXUS_TREE -s COMPLETE ALN -lh HUMAN_LOCATION_LABEL -ph SAMPLNG_PROB_HIGH -pl SAMPLING_PROB_LOW
    # Run subsample_trees.py which is ready now
# We need a script to estimate trees
    # Run find_monophyletic.py to calculate importation events  


#for i in $replicaes; do
#    sed -i.bak 's/OUTPUT_FILE_NAME/$i/g' MPXV_speciesjump_common.xml
 #   done