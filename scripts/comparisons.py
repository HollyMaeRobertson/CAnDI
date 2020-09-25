"""This module contains functions for comparing biparts to each other."""

from objects import Bipart, Rel, Node

import make_trees
import read_trees
import sys

#In an array of arrays finds the smallest array
#that has all elements of an array in question
def get_smallest_match(bipart_array,ingroup):
	
	sort_array = []
	for bipart in bipart_array:
		if(set(ingroup).issubset(set(bipart))):
			sort_array.append(bipart)
	sort_array.sort(key=len)
	return sort_array[0]

def unique_array(tree1, tree2):
    """This function takes three lists of species names ('array' should be
    the sum of tree1 and tree2. It checks to see which species are in one 
    tree but not the other. Returns two lists of the species missing from 
    tree1 (mis1) and the species missing from tree2 (mis2)."""

    # Make a set of all the species in both trees, so we only have unique
    # species.
    array = tree1 + tree2
    all_species = set()
    for x in array:
        all_species.add(x)

    # Self-explanatory.
    mis1 = list(set(all_species) - set(tree1))
    mis2 = list(set(all_species) - set(tree2))

    return mis1, mis2


def bipart_relationship(bp1, bp2):
    """This function takes two biparts (bipart_proper attribute of Bipart()
    class) and compares them. It returns a string describing their 
    relationship"""

    # The species that are in one not the other tell us about the relationship.
    diff1 = list(set(bp1) - set(bp2))
    diff2 = list(set(bp2) - set(bp1))

    # If either bipart refers to a tip (i.e. contains only one species).
    if len(bp1) == 1 or len(bp2) == 1:
        return "uninformative"

    # If they have no taxa at all in common, they're not comparable.
    elif len(diff1) == len(bp1):
        return "no_comp"

    # Concordant means they are completely identical.
    elif len(bp1) == len(bp2) and len(diff1) == 0:
        return "concordant"

    # We shouldn't compare nested biparts because there'll be another bipart
    # somewhere else on the tree which allows a direct comparison.
    elif len(diff1) == 0:
        return "1 nested in 2"
    elif len(diff2) == 0:
        return "2 nested in 1"

    # The only other possibility is conflict.
    else:
        return "conflict"


def comp_biparts(tree1, tree2, name_array1, name_array2, log_name, cutoff, mode, gene_name):
    """This function takes two lists of biparts (tree1 and tree2), the 
    corresponding lists of names (name_array 1 and 2), the name of a log to
    write out to and the bootstrap cutoff value (we ignore nodes below the 
    cutoff)	and returns a list of all concordant or conflicting instances of
    the Rel() class, which contains the unique identifiers of two nodes and 
    the relationship between them.

    NB: The variables tree1 and tree2 are NOT equivalent. tree1 should be 
    the tree you want to map your relationships back onto. Similarly,
    name_array1 should correspond with tree1 and name_array2 with tree2.
    """


    # Initialise things.
    relationship_list = []
    test_bp1 = []
    test_bp2 = []
    rel = ""
    outf = open(log_name + ".log", "w")
    count = 0

    # We should get missing taxa to allow us to exclude them later
    if mode != "s":
        mis1, mis2 = unique_array(name_array1, name_array2)
    else:
        mis1 = []
        mis2 = []

    # The 'tree to map back onto' bipartitions.
    for bp1 in tree1:

        various_relationships = []

        # Make the first bipart to test, removing missing taxa.
        test_bp1 = list(set(bp1.bipart_proper) - set(mis2))
        bp1_len = len(list(bp1.bipart_proper))

        # Removing missing taxa can generate a lot of duplicate
        # bipartitions. We want to only use the most "downstream" of
        # these, which should be the one with the shortest overall
        # bipart proper.
        shortest_bp = True
        for bp in tree1:
            bp_len = len(bp.bipart_proper)
            test_bp = list(set(bp.bipart_proper) - set(mis2))

            if sorted(test_bp) == sorted(test_bp1) and bp_len < bp1_len:
                shortest_bp = False

        if not shortest_bp:
            count += 1
            continue

        for bp2 in tree2:
            # Removing missing taxa. This should only be important
            # in "r" mode.
            test_bp2 = list(set(bp2.bipart_proper) - set(mis1))

            # Assuming the label is a bootstrap/confidence value
            # (out of 100), here we treat it as a cutoff to decide
            # whether or not to include these nodes in the analysis.
            # If the label isn't an integer, we always include it.
        
            #if str.isdigit(bp1.label):
            #    cutoff1 = int(bp1.label)
            #else:
           	    #cutoff1 = 100

            if str.isdigit(bp2.label):
                cutoff2 = int(bp2.label)
            else:
                cutoff2 = 100

            #Only compare the gene biparts, reverse mode does something a bit different
            #so avoid that
            if mode == "n" and cutoff2 < cutoff and cutoff != 0:
                pass

            else:
                # The relationship between the biparts.
                rel = bipart_relationship(test_bp1, test_bp2)
                outf.write(str(rel) + ": " + str(bp1.bipart_proper) +
                           " | " + str(bp2.bipart_proper) + "\n")
                # We only record these two cases.
                if rel == "conflict" or rel == "concordant":
                    if mode == "n" or mode == "s":

                        relation = Rel(rel, bp1.unique_id, bp2.unique_id, gene_name)
                        relation.add_species_bipart(bp1.bipart_proper)
                        relation.add_ortholog_bipart(bp2.bipart_proper)
                    elif mode == "r":
                        relation = Rel(rel, bp2.unique_id, bp1.unique_id)
                        relation.add_species_bipart(bp2.bipart_proper)
                        relation.add_ortholog_bipart(bp1.bipart_proper)
                    various_relationships.append(relation)

                # In a rare subset of cases, directly up from a
                # conflict there is an identical species bipart
                # when missing taxa are removed.
                # This is also in conflict and should be
                # recorded as such.
                if rel == "conflict" and mode == "n":
                    for bp in tree1:
                        test_bp = list(set(bp.bipart_proper) - set(mis2))
                        if sorted(test_bp) == sorted(test_bp1):
                            new_rel = Rel(
                                "conflict", bp.unique_id, bp2.unique_id, gene_name)
                            new_rel.add_species_bipart(bp.bipart_proper)
                            new_rel.add_ortholog_bipart(bp2.bipart_proper)
                            various_relationships.append(new_rel)

        relationship_list.extend(various_relationships)
        count += 1

    return relationship_list


def compare_trees(species_biparts, species_name_array, subtree, extra_names, family, mode, log_name, cutoff, gene_name):
    """This function compares a subtree to a species tree which has already	
    been split into biparts.
    """

    conflicts = []
    concordances = []
    keepgoing = True
    current_node = subtree
    new_names = []

    # We need to add the names that are 'next-door' to this subtree on the
    # overall tree. (Why??? Come back to this when troubleshooting.)

    while True:
        parent = current_node.parent
        if parent != None:
            make_trees.label_duplications(parent, recursive=False)
            if parent.label == "D":
                current_node = parent
            else:
                for i in parent.children:
                    if i != current_node:
                        new_names = read_trees.postorder3(i)
                        new_names = new_names.bipart_proper
                break
        else:
            break

    gene_name_array = read_trees.postorder3(subtree)
    gene_name_array = gene_name_array.bipart_proper
    gene_name_array.extend(new_names)
    
    #If this is a homolog then the earliest node to not have a duplication must be found
    if family == "homolog":
    	gene_name_array.extend(extra_names)

    		

    # Actually making the comparisons.
    subtree_biparts = read_trees.postorder2(subtree, subtrees=True)

    if mode == "n" or mode == "s":
        rels = comp_biparts(species_biparts, subtree_biparts,
                            species_name_array, gene_name_array, log_name, cutoff, mode, gene_name)
    elif mode == "r":
        rels = comp_biparts(subtree_biparts, species_biparts,
                            gene_name_array, species_name_array, log_name, cutoff, mode, gene_name)

    for rel in rels:
        if rel.relation == "conflict":
            conflicts.append(rel)
        elif rel.relation == "concordant":
            concordances.append(rel)

    return conflicts, concordances


def filter_conflicts(conflicts):
    """The function comp_biparts returns a list of *all* conflicts a subtree
    has with the species tree. This includes a large number of redundant 
    conflicts as unlike with concordance (which will always only happen once
    per speices node per gene tree) a gene tree can conflict with multiple 
    nodes on the species tree and multiple times with each node. This 
    function picks one conflict per species node (where one exists), always
    choosing the most upstream place that conflict has ever been recorded on
    both the species and gene trees."""

    # We want to end up with exactly 1 or 0 conflicts for every node present.
    species_nodes = set()
    for conflict in conflicts:

        species_nodes.add(conflict.species_node)
    conflicts_to_return = []

    for node in species_nodes:
        current_conflicts = []

        # We need specifically conflicts that refer to this node.
        for conflict in conflicts:
            if conflict.species_node == node:
                current_conflicts.append(conflict)

        # If there are multiple conflicts, we then take out the ones
        # with the longest ortholog_bipart.
        if len(current_conflicts) == 1:
            correct_conflict = current_conflicts[0]

        else:
            correct_conflict = best_conflict_machine(current_conflicts, [])
            alt_conflict = best_conflict_machine(
                current_conflicts, [correct_conflict])
            if alt_conflict:
                correct_conflict.add_alt_conflict(alt_conflict.ortholog_bipart)
        conflicts_to_return.append(correct_conflict)

    return conflicts_to_return


def best_conflict_machine(current_conflicts, best_conflicts):
    """The 'best' conflict from a list of potential conflicts has the 
    longest ortholog bipart and isn't nested in any other conflicts.
    """
    best_conflict = None
    length = 0
    for conflict in current_conflicts:
        if len(conflict.ortholog_bipart) >= length:
            include_conflict = True
            for conflict2 in best_conflicts:
                rel = bipart_relationship(
                    conflict.ortholog_bipart, conflict2.ortholog_bipart)
                if rel == "1 nested in 2" \
                        or rel == "2 nested in 1" \
                        or rel == "concordant":
                    include_conflict = False
            if include_conflict:
                best_conflict = conflict
                length = len(conflict.ortholog_bipart)
    return best_conflict


def filter_conflicts_for_csv(conflicts):

    # First we're gonna get a list of all the unique species nodes we've got
    # ordered "biggest" to "smallest" in terms of length of bipart.
    species_nodes = set()
    for conflict in conflicts:
        quick_string = str(conflict.species_node) + ":" + \
            str(len(conflict.species_bipart))
        species_nodes.add(quick_string)

    species_nodes = list(species_nodes)
    species_nodes_new = []
    for node in species_nodes:
        correct_node = node.split(":")
        species_nodes_new.append(correct_node)
    species_nodes = species_nodes_new

    def second_thing(list_two_things):
        return int(list_two_things[1])

    species_nodes = sorted(species_nodes, key=second_thing, reverse=True)
    best_conflicts = []

    for node in species_nodes:
        node_actual = node[0]
        current_conflicts = []

        # We only look at conflicts which apply to this species node.
        for conflict in conflicts:
            if conflict.species_node == node_actual:
                current_conflicts.append(conflict)

        # We want to add only the longest conflicts that aren't nested
        # in any others.
        while True:
            best_conflict = best_conflict_machine(
                current_conflicts, best_conflicts)
            if best_conflict:
                best_conflicts.append(best_conflict)
            else:
                break
    return best_conflicts
