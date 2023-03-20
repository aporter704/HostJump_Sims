# %%
import dendropy as dp
import numpy as np
import re, sys, os
import pandas as pd
import subprocess, sqlite3

# %%
def annotate_tip_labels_multitype(tree):
    """Annotate the tips of a multitype tree according to the _location_ label from 
    MASTER

    Args:
        tree (dendropy.Tree): This should be a multitype tree from MASTER 

    Returns:
        list: Element 0 is the tree with annotated tips and element 1 is a 
        pandas data frame with the old and new name.
    """
    taxon_translation_list = []
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            location = node.annotations.get_value('location')
            new_name = node.taxon.label + '_type' + location + '_' + str(round(node.distance_from_root(), 2))
            taxon_translation_list.append([node.taxon.label, new_name])
            node.taxon.label = new_name

    taxon_translation_list = pd.DataFrame(taxon_translation_list)
    taxon_translation_list.columns = ['old_name', 'new_name']
    
    return [tree, taxon_translation_list]

# %%
def annotate_tip_labels_from_table(tree, taxon_translation_list):
    """Annotate the tips of any tree by matching them to a pandas data frame
        with two colums, old_name and new_name, that can be generated with annotate_tip_labels_multitype

    Args:
        tree dendropy.Tree: tree from dendropy
        taxon_translation_list pandas.DataFrame: data frame with two columns, 
            one for the current taxon label and one for the new one to be assigned

    Returns:
        dendropy.Tree : tree with new tip labels
    """
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            temp =  taxon_translation_list.loc[taxon_translation_list.old_name == node.taxon.label, 'new_name']
            node.taxon.label = temp.values[0]
            
    return tree

# %%
"""
multitype_tree = dp.Tree.get_from_path("OUTPUT_FILE_NAME.nexus.tree", schema = "nexus")
annotated_tree_data = annotate_tip_labels_multitype(multitype_tree)
"""

# %%
"""
plain_tree = dp.Tree.get_from_path("OUTPUT_FILE_NAME.newick.tree", schema = "newick")
annotated_plain_tree = annotate_tip_labels_from_table(plain_tree, annotated_tree_data[1])
annotated_plain_tree
"""


# %%
def count_migrations_multitype_tree(tree):
    """Count number of migration events for a multitype tree

    Args:
        tree (dendropy.Tree): tree with internal node annotations. 
            The _migration_ reaction must be stored in node annotations.

    Returns:
        pandas.DataFrame: a data frame with 4 columns, for the reaction type 
            (always migration), the ancestor type, the child type, and the age
            of the event.
    """
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
    
    return migration_events

# %%
"""
count_migrations_multitype_tree(multitype_tree).head()
"""

# %%
def get_tip_ages(tree):
    """Get the age for tips. The age is the distance from the root, 
        and not to be confused with the height.

    Returns:
        pandas.DataFrame: data frame with 3 columns, the tip_label, type, and date (age).
    """
    tip_ages = []
    for i in tree.leaf_node_iter():
        tip_ages.append([i.taxon.label, i.annotations.get_value('location'), i.distance_from_root()])

    tip_ages = pd.DataFrame(tip_ages)
    tip_ages.columns = ['tip_label', 'type', 'date']
    
    return tip_ages

# %%
"""
get_tip_ages(multitype_tree).head()
"""

# %%
def select_tips_to_remove(tree, location_to_retain, location_to_remove, sampling_proportion, opportunistic = False, sampling_start_age = 0):
    """Sample tips according to their _location_ annotation (here location_to_retain)

    Args:
        tree (dendropy.Tree): A multitype tree with _location_ annotations.
        location_to_retain (int): the label to match the _location_ annotation.
        location_to_remove (int): the label to match the _location_ annotation.
        sampling_proportion (float): a number between 0 and 1. The sampling 
            probability of any samples not belonging to location_to_retain
        opportunistic (bool): whether to do opportunistic sampling 
            (after a certain time) nor sample homogeneously over time.
        sampling_start_age (float): the age (from root) after which to start 
            sampling.

    Returns:
        list: A list of tips to be removed.
    """
    location_to_retain = str(location_to_retain)
    if location_to_remove is not None:
        location_to_remove = str(location_to_remove)
        
    tips_to_remove = []
    tip_ages = get_tip_ages(tree)
    if(opportunistic and location_to_remove is None):
        for i in range(tip_ages.shape[0]):
            if tip_ages.date[i] < sampling_start_age:
                tips_to_remove.append(tip_ages.tip_label[i])
            elif tip_ages.type[i] != location_to_retain and tip_ages.date[i] >= sampling_start_age:
                if np.random.binomial(1, sampling_proportion) == 0:
                    tips_to_remove.append(tip_ages.tip_label[i])
    elif(opportunistic and location_to_remove is not None):
        for i in range(tip_ages.shape[0]):
            if tip_ages.date[i] < sampling_start_age:
                tips_to_remove.append(tip_ages.tip_label[i])
            elif tip_ages.type[i] == location_to_remove and tip_ages.date[i] >= sampling_start_age:
                if np.random.binomial(1, sampling_proportion) == 0:
                    tips_to_remove.append(tip_ages.tip_label[i])


    else:
        for i in tree.leaf_node_iter():
            location_temp =  i.annotations.get_value('location')
            if location_temp != location_to_retain:
                if np.random.binomial(1, sampling_proportion) == 0:
                    tips_to_remove.append(i.taxon.label)
                
    return tips_to_remove                
            


# %%
"""
tips_to_remove = select_tips_to_remove(multitype_tree, location_to_retain=1, sampling_proportion= 0.9)
tree_sampled = multitype_tree.clone(depth = 1)
tree_sampled.prune_taxa_with_labels(tips_to_remove) # note that the output tree here is not a mutitype!
tree_sampled.write_to_path("subsampled.tree", "nexus")
"""

# %%
"""
tip_ages = get_tip_ages(multitype_tree)
first_type1 = tip_ages.date[tip_ages.type == "1"].min()
first_type1
"""


# %%
"""
select_tips_to_remove(multitype_tree, location_to_retain=1, sampling_proportion=0.9, opportunistic=True, sampling_start_age=first_type1)
"""


def count_migrations_nonmultitype(tree, target_type):
    """Count the number of putative migration events in a non multitype tree. 
        The number of events are estimated as the number of monophyletic groups
        of a target_type.

    Args:
        tree (dendropy.Tree): A fully bifurcating tree and with tip labels 
            that contain _target_type_. Check regexp below as needed
        target_type (str): stirng denoting types to be matched with regexp.

    Returns:
        list: a list where the first item is the number of importation/migration events. 
            The second element is the age of the first importation event, and the last is 
            the age of the first imported sample.
    """
    for i in tree.leaf_node_iter():
        i.taxon.label = re.sub('^_|^ |_$| $', '', i.taxon.label)

    # A verification for when the prunning is very extensive and no importation events can be identifified
    all_tip_labels = [tip.taxon.label for tip in tree.leaf_node_iter()]
    count_nontarget_types = [tip for tip in all_tip_labels if not target_type in tip]
    if len(count_nontarget_types) == 0:
        return [0, 0, 0]
    if len(all_tip_labels) == len(count_nontarget_types):
        return [0, 0, 0]

    visited_nodes = []
    monophyletic_nodes = []
    clades = []

    for int_node in tree.preorder_node_iter():
        if int_node.is_leaf() or int_node in visited_nodes:
            continue
        child_tips = [i.taxon.label for i in int_node.leaf_iter()]
        match_target = [target_type in i for i in child_tips]
        if all(match_target):
            monophyletic_nodes.append(int_node)
            clades.append(child_tips)
            for sub_node in int_node.preorder_internal_node_iter():
                visited_nodes.append(sub_node)

    unlisted_clades = []
    for clade in clades:
        for taxon in clade:
            unlisted_clades.append(taxon)

    singletons = [tip for tip in tree.leaf_node_iter() if (not tip.taxon.label in unlisted_clades) and (target_type in tip.taxon.label)]
    node_ages = [i.distance_from_root() for i in monophyletic_nodes]
    tip_branching_ages = [i.distance_from_root()-i.edge_length for i in singletons]
    matching_tip_ages = [i.distance_from_root() for i in tree.leaf_node_iter() if target_type in i.taxon.label]
    if(len(tip_branching_ages) == 0):
        tip_branching_ages = 0
    if(len(node_ages) == 0):
        node_ages = 0
    first_import_age = np.min((np.min(tip_branching_ages), np.min(node_ages)))

    return([len(clades)+len(singletons), first_import_age, np.min(matching_tip_ages)])

# %%
def count_migrations_nonmultitype_deprecated(tree, target_type):
    """Count the number of putative migration events in a non multitype tree. 
        The number of events are estimated as the number of monophyletic groups
        of a target_type.

    Args:
        tree (dendropy.Tree): A fully bifurcating tree and with tip labels 
            that contain _target_type_. Check regexp below as needed
        target_type (str): stirng denoting types to be matched with regexp.

    Returns:
        list: a list where the first item is the number of importation/migration events. 
            The second element is the age of the first importation event, and the last is 
            the age of the first imported sample.
    """
    for i in tree.leaf_node_iter():
        i.taxon.label = re.sub('^_|^ |_$| $', '', i.taxon.label)

    # A verification for when the prunning is very extensive and no importation events can be identifified
    all_tip_labels = [tip.taxon.label for tip in tree.leaf_node_iter()]
    count_nontarget_types = [tip for tip in all_tip_labels if not target_type in tip]
    if len(count_nontarget_types) == 0:
        return [0, 0, 0]
    if len(all_tip_labels) == len(count_nontarget_types):
        return [0, 0, 0]
    if len(count_nontarget_types) < 10:
        return [0, 0, 0]


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

    return [len(importation_nodes), np.min(importation_ages), np.min(tip_ages)]


# %%
"""
print( count_migrations_nonmultitype(plain_tree, "type1") )
out = get_tip_ages(plain_tree)
is_type1 = [i.taxon.label for i in plain_tree.leaf_node_iter() if len(re.findall("type1", i.taxon.label)) > 0]
is_type1 = [i in is_type1 for i in out.tip_label]
out.loc[is_type1, ].tail(10)
"""

# %%
"""
out = count_migrations_multitype_tree(multitype_tree)
"""

# %%
"""
cond1 = out.reaction == "Migration" 
cond2 = out.child_type == "1"
len([i for i in range(len(cond1)) if cond1[i] and cond2[i]])
"""

# %%
def simulate_tree_master(path_to_beast2, path_to_template_file, prefix):
    """ Call MASTER to simulate trees from within python using a template xml and replacing 
        OUTPUT_FILE_NAME by a prefix to idenitfy easily.

    Args:   
        path_to_beast2 (str)): a string with the full path to BEAST2 with MASTER installed.
        path_to_template_file (str): a string with the template of the file with the template
        prefix (str): a prefix for the output file name
    """
    template_xml = open(path_to_template_file, 'r').read()
    replicate_xml = re.sub('OUTPUT_FILE_NAME', prefix, template_xml)
    open(prefix + '.xml', 'w').writelines(replicate_xml)
    os.system(path_to_beast2 + ' ' +prefix+'.xml')

# %%
"""
path_to_beast2 = "/home/sebastiand/phyloApps/beast2.6.3/bin/beast"
path_to_template_file = "/home/sebastiand/Dropbox/projects_WORKING/speciesJumpSimulations/HostJump_Sims/master_scripts/final_EBOV_speciesjump_rare.xml"
prefix = "testing_function"
simulate_tree_master(path_to_beast2 = path_to_beast2, path_to_template_file = path_to_template_file, prefix = prefix)
"""

# %%


# %%
# We'll also need some code to push results to a sqlite data base. This may not  need to be a function afterall


