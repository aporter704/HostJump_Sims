#!/usr/bin/env python3

import os, sys, re

replicates = [i for i in range(5)] # number of replicates to run

template = 'MPXV_speciesjump_common.xml' # Template xml
prefix = 'set1_MPXV_common_r' # the prefix we want to use


for i in replicates[:3]:
    replicate_xml = open(template, 'r').read()
    output_name = prefix+str(i)
    replicate_xml = re.sub('OUTPUT_FILE_NAME', output_name, replicate_xml)
    open(output_name + '.xml', 'w').writelines(replicate_xml)
    os.system('~/phyloApps/beast2.6.3/bin/beast ' + output_name+'.xml')
    os.system('./treePostProcessing.py -t ' + output_name + '.nexus.tree -lh 1')
# Here execute code to simulate the sequence alignment   



#for i in $replicaes; do
#    sed -i.bak 's/OUTPUT_FILE_NAME/$i/g' MPXV_speciesjump_common.xml
 #   done