#!/usr/bin/env python3

# Requires subsample_trees.py, find_monophyletic.py, tree_post_processing.py, make_lsd_dates.py

import os, sys, re, subprocess, sqlite3

# The parameters below will set ups simulations:
template = 'test_template.xml' # template for running simulation, usually a MASTER file 
prefix = 'set_test_template_r' # the prefix for output files
beast_path = '~/phyloApps/beast2.6.3/bin/beast' # Path to beast
seqgen_path = '~/phyloApps/Seq-Gen-1.3.4/source/seq-gen' # path to seq-gen
clock_rate = '1e-4' # clock rate in subs/site/unit time
genome_size = '30000' # genome size in nt
iqtree_path = '~/phyloApps/iqtree-2.1.3-Linux/bin/iqtree2' # path to iqtree 
lsd_path = '~/phyloApps/lsd-0.3beta/bin/lsd_unix' # path to lsd
target_type = 'type1' # type label to match. In tree_post_processing.py we match the 'location' tag
num_replicates = 3 # the number of replicates to run
database_name = 'test_template_results.db'

for i in range(num_replicates):
    replicate_xml = open(template, 'r').read()
    output_name = prefix+str(i)
    replicate_xml = re.sub('OUTPUT_FILE_NAME', output_name, replicate_xml)
    open(output_name + '.xml', 'w').writelines(replicate_xml)
    os.system(beast_path + ' ' +output_name+'.xml')
    tree_processing = re.split(' |\n', os.popen('./tree_post_processing.py -t ' + output_name + '.nexus.tree -lh 1 -rc '+clock_rate+' -gl '+genome_size+' -seqGen ' + seqgen_path).read()) # manipulate trees

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
        estimated_trees.append(aln+'.treefile.result.date.nexus')

# Make table and compile statistics:
    con = sqlite3.connect(database_name)
    cur = con.cursor() 
    table_name = re.sub('_nMig.+', '', sim_tree_file)
    table_name
    cur.execute('CREATE TABLE IF NOT EXISTS '+table_name+'(tree_name, tmrca, jumps, first_migration_age, first_detected_case)')
    con.commit()
    
    # Simulated painted tree *
    tmrca = re.sub('[.]nexus.+', '', re.sub('.+splitAge', '', sim_tree_file))
    num_migrations = re.sub('_.+', '', re.sub('.+_nMig', '', sim_tree_file))
    first_migration_age = re.sub('_.+', '', re.sub('.+migAge', '', sim_tree_file))
    first_detected_case = '\'NA\''
    summary_stats = tmrca+','+num_migrations+','+first_migration_age+','+first_detected_case

    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'simulated_tree\','+summary_stats+')')
    con.commit()    

    # Simulated non painted tree *
    non_painted_tree = re.sub('nexus', 'newick', sim_tree_file)
    sim_non_painted_tree_stats = os.popen('./find_monophyletic.py -t '+non_painted_tree+' -f newick -l '+target_type).readlines()
    sim_non_painted_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(sim_non_painted_tree_stats) if i[0] in [1, 3, 5, 7]]
    sim_non_painted_tree_stats = ','.join(sim_non_painted_tree_stats)
    sim_non_painted_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'simulated_non_painted_tree\','+sim_non_painted_tree_stats+')')
    con.commit()
      
    # Simulated pruned high sampling *
    high_sampling_tree = re.sub('[.]nexus.+', '_prunedHighSamp.nexus.tree', sim_tree_file)
    high_sampling_tree_stats = os.popen('./find_monophyletic.py -t '+high_sampling_tree+' -f nexus -l '+target_type).readlines()
    high_sampling_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(high_sampling_tree_stats) if i[0] in [1, 3, 5, 7]]
    high_sampling_tree_stats = ','.join(high_sampling_tree_stats)
    high_sampling_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'simulated_high_sampling_tree\','+high_sampling_tree_stats+')')
    con.commit()

    # Simulated pruned low sampling *
    low_sampling_tree = re.sub('[.]nexus.+', '_prunedLowSamp.nexus.tree', sim_tree_file)
    low_sampling_tree_stats = os.popen('./find_monophyletic.py -t '+low_sampling_tree+' -f nexus -l '+target_type).readlines()
    low_sampling_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(low_sampling_tree_stats) if i[0] in [1, 3, 5, 7]]
    low_sampling_tree_stats = ','.join(low_sampling_tree_stats)
    low_sampling_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'simulated_low_sampling_tree\','+low_sampling_tree_stats+')')
    con.commit()
    
    # Simulated opportunistic sampling *
    opp_sampling_tree = re.sub('[.]nexus.+', '_prunedOppSamp.nexus.tree', sim_tree_file)
    opp_sampling_tree_stats = os.popen('./find_monophyletic.py -t '+opp_sampling_tree+' -f nexus -l '+target_type).readlines()
    opp_sampling_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(opp_sampling_tree_stats) if i[0] in [1, 3, 5, 7]]
    opp_sampling_tree_stats = ','.join(opp_sampling_tree_stats)
    opp_sampling_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'simulated_opp_sampling_tree\','+opp_sampling_tree_stats+')')
    con.commit()

    # Estimated complete tree *
    estimated_complete_tree = estimated_trees[0]
    estimated_complete_tree_stats = os.popen('./find_monophyletic.py -t '+estimated_complete_tree+' -f nexus -l '+target_type).readlines()
    estimated_complete_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(estimated_complete_tree_stats) if i[0] in [1, 3, 5, 7]]
    estimated_complete_tree_stats = ','.join(estimated_complete_tree_stats)
    estimated_complete_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'estimated_complete_tree\','+estimated_complete_tree_stats+')')
    con.commit()

    # Estimated pruned high sampling * 
    estimated_high_sampling_tree = estimated_trees[3]
    estimated_high_sampling_tree
    estimated_high_sampling_tree_stats = os.popen('./find_monophyletic.py -t '+estimated_high_sampling_tree+' -f nexus -l '+target_type).readlines()
    estimated_high_sampling_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(estimated_high_sampling_tree_stats) if i[0] in [1, 3, 5, 7]]
    estimated_high_sampling_tree_stats = ','.join(estimated_high_sampling_tree_stats)
    estimated_high_sampling_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'estimated_high_sampling_tree\','+estimated_high_sampling_tree_stats+')')
    con.commit()

    # Estimated pruned low sampling *
    estimated_low_sampling_tree = estimated_trees[2]
    estimated_low_sampling_tree
    estimated_low_sampling_tree_stats = os.popen('./find_monophyletic.py -t '+estimated_low_sampling_tree+' -f nexus -l '+target_type).readlines()
    estimated_low_sampling_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(estimated_low_sampling_tree_stats) if i[0] in [1, 3, 5, 7]]
    estimated_low_sampling_tree_stats = ','.join(estimated_low_sampling_tree_stats)
    estimated_low_sampling_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'estimated_low_sampling_tree\','+estimated_low_sampling_tree_stats+')')
    con.commit()

    # Estimated opportunistic sampling
    estimated_opp_sampling_tree = estimated_trees[1]
    estimated_opp_sampling_tree
    estimated_opp_sampling_tree_stats = os.popen('./find_monophyletic.py -t '+estimated_opp_sampling_tree+' -f nexus -l '+target_type).readlines()
    estimated_opp_sampling_tree_stats = [re.sub('\n', '', i[1]) for i in enumerate(estimated_opp_sampling_tree_stats) if i[0] in [1, 3, 5, 7]]
    estimated_opp_sampling_tree_stats = ','.join(estimated_opp_sampling_tree_stats)
    estimated_opp_sampling_tree_stats
    cur.execute('INSERT INTO '+table_name+' VALUES ('+'\'estimated_opp_sampling_tree\','+estimated_opp_sampling_tree_stats+')')
    con.commit()

