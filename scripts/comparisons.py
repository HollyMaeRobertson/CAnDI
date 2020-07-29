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
		lengths = []
		
		# Make the first bipart to test, removing missing taxa
		test_bp1 = list(set(bp1.bipart_proper) - set(mis2))

		# Removing missing taxa can generate a lot of duplicate 
		# bipartitions. We want to only use the most "downstream" of 
		# these, which should be the last one in the list, so if there
		# is an identical bipart further on in the list tree1, we pass 
		# the others until we reach that one.
		remaining_bps = [set(i.bipart_proper) - set(mis2) for i in tree1[count+1:]]
		if test_bp1 in remaining_bps:
			count += 1
			continue

		for bp2 in tree2:
			
			# Removing missing taxa as above.
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
					elif mode == 'r':
						relation = Rel(rel, bp2.unique_id, bp1.unique_id)
					relation.add_species_bipart(test_bp1)
					relation.add_ortholog_bipart(test_bp2)
					various_relationships.append(relation)
					lengths.append(len(test_bp2))

		# If there was only one relation in the subtree that was 
		# 'conflict' or 'concordant.
		if len(various_relationships) == 1:
			relation = various_relationships[0]
			relationship_list.append(relation)
		
		# If more than one node in the subtree was 'conflict' or 
		# 'concordant', we count all the nodes that were longest???
		elif len(various_relationships) != 0:
			num = 1000000000000000
			counter = 0
			
			# Why did I do this??? Come back to this. 
			for k in lengths:
				if k < num:
					num = k
					relation = various_relationships[counter]
				counter += 1
				relationship_list.append(relation)

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



