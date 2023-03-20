from get_tree_stats import *
from math import *

path_to_beast2 = "/home/sebastiand/phyloApps/beast2.6.3/bin/beast"
path_to_template_file = "/home/sebastiand/Dropbox/projects_WORKING/speciesJumpSimulations/HostJump_Sims/master_scripts/final_HeV_speciesjump_common.xml"
database_name = "HeV.db"
dataset = "HeV_common"

for i in range(1500):
    prefix = dataset+str(i)
    simulate_tree_master(path_to_beast2 = path_to_beast2, path_to_template_file = path_to_template_file, prefix = prefix)

    test_lines = open(prefix+".nexus.tree").read()
    if len(test_lines) < 50:
        continue

    multitype_tree = dp.Tree.get_from_path(prefix+".nexus.tree", schema = "nexus")

    if len(multitype_tree.edges()) < 10:
        continue

    is_migration = [node.annotations.get_value("reaction") == "Migration" for node in multitype_tree.postorder_node_iter()]
    if not any(is_migration):
        continue

    tip_labels_type = annotate_tip_labels_multitype(multitype_tree)
    multitype_tree = tip_labels_type[0]

    nontyped_tree = dp.Tree.get_from_path(prefix+".newick.tree", schema = "newick")
    nontyped_tree = annotate_tip_labels_from_table(nontyped_tree, tip_labels_type[1])
    nontyped_tree.write_to_path(prefix+"_nontyped.newick.tree", "newick")


    tip_labels = [tip.taxon.label for tip in multitype_tree.leaf_node_iter()]
    tips_with_label = [matching_tip for matching_tip in tip_labels if len(re.findall("type0", matching_tip))]
    if len(tips_with_label) < 4:
        continue

    # Prunning tree with high, low, opportunistic0.4, and opportunistic0.05
    prune_high_prop = select_tips_to_remove(multitype_tree, "2", 0.5, opportunistic=False)
    tree_prune_high_prop = nontyped_tree.clone(depth = 1)
    tree_prune_high_prop.prune_taxa_with_labels(prune_high_prop) 
    tree_prune_high_prop.write_to_path(prefix+"_high_prop.newick.tree", "newick")

    prune_low_prop = select_tips_to_remove(multitype_tree, "2", 0.05, opportunistic=False)

    if len(prune_low_prop) >= (len(tip_labels)-5):
        continue

    tree_prune_low_prop = nontyped_tree.clone(depth = 1)
    tree_prune_low_prop.prune_taxa_with_labels(prune_low_prop)    
    tree_prune_low_prop.write_to_path(prefix+"_low_prop.newick.tree", "newick")

    tip_ages = get_tip_ages(multitype_tree)
    first_target_date = tip_ages.date[tip_ages.type == "2"].min()

    prune_opportunistic_high = select_tips_to_remove(multitype_tree, location_to_retain = "2", sampling_proportion = 0.5, opportunistic=True, sampling_start_age = first_target_date)
    tree_prune_opportunistic_high = nontyped_tree.clone(depth = 1)
    tree_prune_opportunistic_high.prune_taxa_with_labels(prune_opportunistic_high)
    tree_prune_opportunistic_high.write_to_path(prefix+"_opp_high.newick.tree", "newick")

    prune_opportunistic_low = select_tips_to_remove(multitype_tree, location_to_retain = "2", sampling_proportion = 0.05, opportunistic=True, sampling_start_age = first_target_date)
    tree_prune_opportunistic_low = nontyped_tree.clone(depth = 1)
    tree_prune_opportunistic_low.prune_taxa_with_labels(prune_opportunistic_low)
    tree_prune_opportunistic_low.write_to_path(prefix+"_opp_low.newick.tree", "newick")

    true_num_migration_events = count_migrations_multitype_tree(multitype_tree)
    from0_to_1 = [true_num_migration_events.ancestor_type[i] == "0" and true_num_migration_events.child_type[i] == "1" for i in range(true_num_migration_events.shape[0])]
    true_migration_0_1 = true_num_migration_events.loc[from0_to_1, ].shape[0]

    from0_to_2 = [true_num_migration_events.ancestor_type[i] == "0" and true_num_migration_events.child_type[i] == "2" for i in range(true_num_migration_events.shape[0])]
    true_migration_0_2 = true_num_migration_events.loc[from0_to_2, ].shape[0]
    from1_to_2 = [true_num_migration_events.ancestor_type[i] == "1" and true_num_migration_events.child_type[i] == "2" for i in range(true_num_migration_events.shape[0])]
    true_migration_1_2 = true_num_migration_events.loc[from1_to_2, ].shape[0]
    true_first_migration_age0 = true_num_migration_events.time[from0_to_2].min()
    if isnan(true_first_migration_age0):
        true_first_migration_age0 = -1
    true_first_migration_age1 = true_num_migration_events.time[from1_to_2].min()
    if isnan(true_first_migration_age1):
        true_first_migration_age1 = -1

    nontype_migration_0_1 = count_migrations_nonmultitype(nontyped_tree, "type2")
    nontype_num_migrations_0_1 = nontype_migration_0_1[0]
    nontype_first_migration_0_1 = nontype_migration_0_1[1]

    prune_high_prop_migration_0_1 = count_migrations_nonmultitype(tree_prune_high_prop, "type2")
    prune_high_num_migrations_0_1 = prune_high_prop_migration_0_1[0]
    prune_high_first_migration_0_1 = prune_high_prop_migration_0_1[1]

    prune_low_prop_migration_0_1 = count_migrations_nonmultitype(tree_prune_low_prop, "type2")
    prune_low_num_migrations_0_1 = prune_low_prop_migration_0_1[0]
    prune_low_first_migration_0_1 = prune_low_prop_migration_0_1[1]

    prune_opportunistic_high_migration_0_1 = count_migrations_nonmultitype(tree_prune_opportunistic_high, "type2")
    prune_opportunistic_high_num_migrations_0_1 = prune_opportunistic_high_migration_0_1[0]
    prune_opportunistic_high_first_migration_0_1 = prune_opportunistic_high_migration_0_1[1]

    prune_opportunistic_low_migration_0_1 = count_migrations_nonmultitype(tree_prune_opportunistic_low, "type2")
    prune_opportunistic_low_num_migrations_0_1 = prune_opportunistic_low_migration_0_1[0]
    prune_opportunistic_low_first_migration_0_1 = prune_opportunistic_low_migration_0_1[1]

    output_data = [] # Here store all the output data and convert nan to 0

    # Prune the trees and get stats
    con = sqlite3.connect(database_name)
    cur = con.cursor() 
    table_name = dataset 
    cur.execute("CREATE TABLE IF NOT EXISTS "+table_name+"(tree_name, true_migration_from0_to1, true_migration_from0_to2, true_migration_from1_to2, true_first_migration_age0, true_first_migration_age1, nontyped_migration_to2, nontyped_first_migration_to2, prune_high_num_migrations, prune_high_first_migration_0_1, prune_low_num_migrations, prune_low_first_migration_0_1, prune_opportunistic_high_num_migrations_0_1, prune_opportunistic_high_first_migration_0_1, prune_opportunistic_low_num_migrations_0_1, prune_opportunistic_low_first_migration_0_1);")
    con.commit()
    cur.execute("INSERT INTO "+table_name+" VALUES ("+"\""+dataset+str(i)+"\","+str(true_migration_0_1)+","+str(true_migration_0_2)+","+str(true_migration_1_2)+","+str(true_first_migration_age0)+","+str(true_first_migration_age1)+","+str(nontype_num_migrations_0_1)+","+str(nontype_first_migration_0_1)+","+str(prune_high_num_migrations_0_1)+","+str(prune_high_first_migration_0_1)+","+str(prune_low_num_migrations_0_1)+","+str(prune_low_first_migration_0_1)+","+str(prune_opportunistic_high_num_migrations_0_1)+","+str(prune_opportunistic_high_first_migration_0_1)+","+str(prune_opportunistic_low_num_migrations_0_1)+","+str(prune_opportunistic_low_first_migration_0_1)+");")
    con.commit()    
    
os.system("tar -cvzf "+dataset+".tar.gz "+dataset+"*nexus*")
os.system("rm "+dataset+"*nexus*")
os.system("rm "+dataset+"*newick*")
os.system("rm "+dataset+"*json*")
os.system("rm "+dataset+"*xml*")
