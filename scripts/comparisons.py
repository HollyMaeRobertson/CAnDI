"""This module contains functions for comparing biparts to each other."""

from objects import Bipart, Rel, Node

import make_trees, read_trees

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

def comp_biparts(tree1, tree2, name_array1, name_array2, log_name, cutoff, mode):
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
	mis1, mis2 = unique_array(name_array1, name_array2)
	
	# The 'tree to map back onto' bipartitions.
	for bp1 in tree1:
		various_relationships = []
		
		# Make the first bipart to test, removing missing taxa
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

                        if test_bp == test_bp1 and bp_len < bp1_len:
                                shortest_bp = False

		if not shortest_bp:
			count += 1
			continue

		for bp2 in tree2:
			
			# Removing missing taxa as above.
                        # (This shouldn't do much as all taxa should be in the 
                        # species tree though?
			test_bp2 = list(set(bp2.bipart_proper) - set(mis1))

			# Assuming the label is a bootstrap/confidence value
			# (out of 100), here we treat it as a cutoff to decide
			# whether or not to include these nodes in the analysis.
			# If the label isn't an integer, we always include it.
			if str.isdigit(bp1.label):
				cutoff1 = int(bp1.label)
			else:
				cutoff1 = 100

			if str.isdigit(bp2.label):
				cutoff2 = int(bp2.label)
			else:
				cutoff2 = 100

			if cutoff1 < cutoff+1 or cutoff2 < cutoff+1:
				pass
			
			else:
				# The relationship between the biparts.
				rel = bipart_relationship(test_bp1, test_bp2)
				outf.write(str(rel) + ": " + str(bp1.bipart_proper) + " | " + str(bp2.bipart_proper) + "\n")
				
				# We only record these two cases.
				if rel == "conflict" or rel == "concordant":
					if mode == 'n' or 's':
						relation = Rel(rel, bp1.unique_id, bp2.unique_id)
                                                relation.add_species_bipart(bp1.bipart_proper)
                                                relation.add_ortholog_bipart(bp2.bipart_proper)
					elif mode == 'r':
						relation = Rel(rel, bp2.unique_id, bp1.unique_id)
                                                relation.add_species_bipart(bp2.bipart_proper)
                                                relation.add_ortholog_bipart(bp1.bipart_proper)
					various_relationships.append(relation)

		relationship_list.extend(various_relationships)
                count += 1

	return relationship_list

def compare_trees(tree1_biparts, name_array1, tree2, mode, log_name, cutoff):
	"""This function compares a subtree to tree1, which has already	been 
	split into biparts.
	"""

	conflicts = []
	concordances = []
	keepgoing = True
	current_node = tree2
	new_names = []
	
	# We need to add the names that are 'next-door' to this subtree on the 
	# overall tree. (Why??? Come back to this when troubleshooting.)
	while True:
		parent = current_node.parent
		if parent != None:	
			make_trees.label_duplications(parent, recursive = False)
			if parent.label == 'D':
				current_node = parent
			else:
				for i in parent.children:
					if i != current_node:
						new_names = read_trees.postorder3(i)
						new_names = new_names.bipart_proper
				break
		else:
			break
		
	name_array2 = read_trees.postorder3(tree2)
	name_array2 = name_array2.bipart_proper
	name_array2.extend(new_names)
	
	# Actually making the comparisons.
	tree2_biparts = read_trees.postorder2(tree2, subtrees = True)
		
	if mode == "n" or "s":
		rels = comp_biparts(tree1_biparts, tree2_biparts, name_array1, name_array2, log_name, cutoff, mode)
	elif mode == "r":
		rels = comp_biparts(tree2_biparts, tree1_biparts, name_array2, name_array1, log_name, cutoff, mode)
	
	for rel in rels:
		if rel.relation == 'conflict':
			conflicts.append(rel)
		elif rel.relation == 'concordant':
			concordances.append(rel)

	return conflicts, concordances


### SHOULD THESE GO IN "analysis.py"???

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
                
                # We want to take out the conflict with only the longest 
                # species_bipart(s).
                length = len(current_conflicts[0].species_bipart)
                filtered_conflicts = []
                for conflict in current_conflicts:
                        if len(conflict.species_bipart) > length:
                                length = len(conflict.species_bipart)
                                filtered_conflicts = []
                                filtered_conflicts.append(conflict)
                        elif len(conflict.species_bipart) == length:
                                filtered_conflicts.append(conflict)
                
                # If there are multiple conflicts with the same length of
                # species_bipart, we then take out the ones with the longest
                # ortholog_bipart
                current_conflicts = filtered_conflicts
                if len(current_conflicts) == 1:
                        correct_conflict = current_conflicts[0]
                else:
                        filtered_conflicts = []
                        length = len(current_conflicts[0].ortholog_bipart)
                        for conflict in current_conflicts:
                                if len(conflict.ortholog_bipart) > length:
                                        length = len(conflict.ortholog_bipart)
                                        filtered_conflicts = []
                                        filtered_conflicts.append(conflict)
                                elif len(conflict.ortholog_bipart) == length:
                                        filtered_conflicts.append(conflict)
                # At this point we just pick a random conflict from those that remain.
                if len(filtered_conflicts) == 1:
                        correct_conflict = filtered_conflicts[0]
                else:
                        correct_conflict = filtered_conflicts[0]
                conflicts_to_return.append(correct_conflict)

        return conflicts_to_return

def best_conflict_machine(current_conflicts, best_conflicts):
        """This function is part of filter_conflicts_for_csv, and it is meant to
        return either the "best" conflict in a list of conflicts, or None, if 
        there is no remaining best conflict.
        """
        best_conflict = None
        length = 0
        for conflict in current_conflicts:
                if len(conflict.ortholog_bipart) >= length:
                        include_conflict = True
                        if best_conflicts:
                                for conflict2 in best_conflicts:
                                        rel = bipart_relationship(conflict, conflict2)
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
                quick_string = str(conflict.species_node) + ":" + str(len(conflict.species_bipart))
                species_nodes.add(quick_string)

        species_nodes = list(species_nodes)
        species_nodes_new = []
        for node in species_nodes:
                correct_node = node.split(":")
                species_nodes_new.append(correct_node)
        species_nodes = species_nodes_new

        def second_thing(list_two_things):
                return int(list_two_things[1])

        species_nodes = sorted(species_nodes, key=second_thing, reverse = True)
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
                        best_conflict = best_conflict_machine(current_conflicts, best_conflicts)
                        if best_conflict:
                                best_conflicts.append(best_conflict)
                        else:
                                break

        return best_conflicts
