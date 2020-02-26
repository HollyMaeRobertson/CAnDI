'''
This program is written to find gene tree conflict and concordance
between a species tree and a set of gene trees (homolog) trees, 
accounting for duplication events. When reverse is True it only takes
one homolog tree.

It can then return either (when reverse is False) two species 
trees (in newick format) with node labels corresponding to the 
number of concordances or  conflicts at that node, 

or (when reverse is True) a single homolog tree with node labels showing 
conflict with the species tree, concordance with the species tree, 
or a duplication event.
'''

import sys, os
from node import Node
 
def build(instr):
	'''This takes in the tree in a newick format, as a string, then 
	puts it in a data structure that can be easily traversed'''

	root = None
	name_array =[]
	index = 0
	nextchar = instr[index]
	begining = "Yep"
	keepgoing = True
	current_node = None
	
	while keepgoing == True:
		# The first node should be a root for our tree
		if nextchar == "(" and begining == "Yep":
				
			root = Node()
			current_node = root
			begining = "No"
		
		# An open bracket in newick format indicates a new node
		elif nextchar == "(" and begining == "No":
		
			newnode = Node()
			current_node.add_child(newnode)
			current_node = newnode
		
		# The commas separate species in the clade so we need to
		# move to the parent node in preparation for the next 
		# related species
		elif nextchar == ',':
		
			current_node = current_node.parent
		
		# Closing brackets close clades, so we need to move up one
		elif nextchar == ")":
			current_node = current_node.parent
			index += 1
			
			# Now we check if there are any labels that apply 
			# to this node written just after the ')'
			
			nextchar = instr[index]
			while True:	
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			
			current_node.label = name
			index -= 1
		
		# The ends of trees in newick format are indicated by semicolons
		# so we stop building the tree
		elif nextchar == ';':
		
			keepgoing = False
			break
		
		# This indicates you have branch lengths so it grabs the branch
		# lengths, turns them into floats and puts them in the current node
		elif nextchar == ":":
			index += 1
			nextchar = instr[index]
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				branch += nextchar
				index += 1
				nextchar = instr[index]
			current_node.length = float(branch)
			index -= 1
		
		# Whitespace means nothing so we just pass it
		elif nextchar == ' ':
			pass
            	
		# Gene trees have locus labels preceded by '@' - we want our tree to 
		# only have taxa labels so we can look for conflict and concordance,
		# but we'll need to put the locus labels back at the end, so we store 
		# them as node attributes
		elif nextchar == '@':
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				
				name += nextchar
				index += 1
				nextchar = instr[index]
			current_node.locus = name
			index -= 1
				
		# External nodes (tips) are different to internal nodes
		else: 
		
			newnode = Node()
			current_node.add_child(newnode)
			current_node = newnode
			current_node.istip = True
			locus = ''
			
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[' or nextchar == '@':
					break

				name += nextchar
				index += 1
				nextchar = instr[index]
					
			current_node.label = name
			name_array.append(name)
			index -= 1
		
		# Moves to the next character
		if index < len(instr) - 1:
			index += 1
		nextchar = instr[index]
		name = ""
		branch = ""
	
	return root,name_array


def postorder3(root, bipart = None):
	'''This traverses the entire tree downstream of the 
	specified root node and puts the names of all the tips
	into the list 'bipart'.'''
	
	# bipart always needs to start empty
	if bipart is None:
		bipart = []
	
	
	for i in root.children:
		if i.istip:
			bipart.append(i.label)
	
		postorder3(i, bipart) 

	return bipart


def postorder2(root, total_list = None, subtrees = False):
	'''This traverses the entire tree downstream of the specified root
	and makes a list (total_list) of all the possible bipartitions, using postorder3.'''
 	
	# total_list should always start off empty
	if total_list is None:
		total_list = []

	# When subtrees is True, we want to make an extra bipartition that encompasss the whole subtree
	if subtrees == True:
		bipart = []
		other_side = [] # There is no other side as this is the whole tree
		
		bipart.append(str(root.length))
		bipart.append(str(root.label))
		bipart.append(other_side)

		bipart = postorder3(root, bipart)
		
		total_list.append(bipart)
	
		subtrees = False

	# Making a bipart to add to total_list
	for i in root.children:
		if i.children:
			bipart = []

			bipart.append(str(i.length)) 
			bipart.append(str(i.label)) 			
				
			# The other side of the bipart is information we want to go into it
			one_up = i.parent
			other_side = []

			for j in one_up.children:
				if j == i:
					pass
				else:
					if j.istip:
						other_side.append(j.label)
					other_side = postorder3(j)
				
			bipart.append(other_side)
		
			# Appends all the tips in the bipart into the list 'bipart'
			bipart = postorder3(i, bipart)
			# Add bipart to the total_list
			total_list.append(bipart)

		postorder2(i, total_list) 		
	
	return total_list


def unique_array(array, tree1, tree2):
	'''Takes three lists of species names and checks to see which species 
	are in one tree but not the other. Returns list of the species missing 
	from tree1 (mis1) and the species missing from tree2 (mis2)'''

	all_species = set()
	
	for x in array:
		all_species.add(x)
	
	mis1 = list(set(all_species) - set(tree1)) 
	mis2 = list(set(all_species) - set(tree2)) 

	return mis1, mis2


def bipart_properties(bp1, bp2):
	'''This function takes two bipart lists made by postorder2 and compares
	them. It returns a string describing their relationship'''
	
	# To make the comparison we need to know their differences
	diff1 = list(set(bp1) - set(bp2))
	diff2 = list(set(bp2) - set(bp1))

	# If the biparts refer to tips (i.e. contain only one species)...
	if len(bp1) == 1 or len(bp2) == 1:
		return "uninformative"

	# If they have nothing at all in common...
	elif len(diff1) == len(bp1):
		return "no_comp"

	# If they are completely identical (concordant)...
	elif len(bp1) == len(bp2) and len(diff1) == 0:
		return "concordant"
	
	# If one is nested in the other...
	elif len(diff1) == 0:
		return "1 nested in 2"
	elif len(diff2) == 0:
		return "2 nested in 1"
	
	# If none of the above cases apply, they must be in conflict.
	else:
		return "conflict"


def make_string(bipart):
	'''This is to make it easier to understand the output when 
	printing things. It makes a neat string out of a bipart 
	made in postorder2 that shows both sides'''
	
	name = ''
	bp = bipart[3:]
	not_bp = bipart[2]

	for i in bp:
		name += i + " "
	
	name += '| '
	
	for i in not_bp:
		name += i + " "

	return str(name)


def comp_biparts(tree1, tree2, name_array1, name_array2, log_name, cutoff):
	'''This function takes two list of biparts (tree1 and tree2), the corresponding 
	lists of names (name_array 1 and 2), a name of a log to write out to, the bootstrap 
	cutoff value (we ignore nodes below the cutoff) and returns a list of all the nodes 
	with concordant or conflicting relationships and which one it was. 
	
	This should be used to compare subtrees (tree2) to the species tree (tree1) most of 
	the time.'''

	# Initialise things
	relationship_list = []
	test_bp1 = []
	test_bp2 = []
	rel = ''
	outf = open(log_name + ".log", "w")
	count = 0

	# Getting missing taxa so we can exclude them
	all_names = name_array1 + name_array2
	mis1,mis2 = unique_array(all_names,name_array1,name_array2)
	
	# Comparing all the biparts pairwise
	for i in tree1:
		various_relationships = []	
		lengths = []

		for j in tree2:
			
			# Remove missing taxa
			test_bp1 = list(set(i[3:]) - set(mis2)) 
			test_bp2 = list(set(j[3:]) - set(mis1))

			# Mostly to make output legible
			not_bp1 = i[2]
			not_bp2 = j[2]
			bp1_string = make_string(i)			
			bp2_string = make_string(j)
			
			# What is the relationship between them?
			rel = bipart_properties(test_bp1, test_bp2)

			# We only want to record the relationship if the bootstrap > cutoff, OR if one or the other node doesn't have a bootstrap
			# Problematic as the label might not actually be a bootstrap even if it is an integer
			if type(i[1]) != int or type(j[1]) != int or str(i[1]) == "" or str(j[1]) == ""  or int(i[1]) > cutoff and int(j[1]) > cutoff:
				
				outf.write(str(rel) + ": " + bp1_string + i[1] + " " + bp2_string + j[1] + "\n")
				
				if rel == "conflict" or rel == "concordant":
					j_info = [str(rel), bp2_string, str(j[1]), str(j[0])]
					various_relationships.append(j_info)
					lengths.append(len(test_bp2))
					
			# Otherwise writes to the log file but doesn't store the relationship in the program
			else:
				outf.write(str(rel) + ": " + str(test_bp1) + i[1] +  " " + str(test_bp2) + j[1] + " "  + "Node not supported" + "\n")
		
		 
		# If there is only one thing in various_relationships, great! 
		# We put it and the node into the relationship_list
		if len(various_relationships) == 1:
			j_info = various_relationships[0]
			print "\n" + str(count) + " " + j_info[0] + "\t" + bp1_string + " " + str(i[1]) + " " + str(i[0]) + "    " + j_info[1] + " " + j_info[2] + " " + j_info[3]
			relationship_list.append([i,j_info[0]])
			
		elif len(various_relationships) != 0:
			num = 100000000000000000
			counter = 0

			for k in lengths:
				if k < num:
					num = k
					j_info = various_relationships[counter]
				counter += 1
			
			print "\n" + str(count) + " " + j_info[0] + "\t" + bp1_string + " " + str(i[1]) + " " + str(i[0]) + "    " + j_info[1] + " " + j_info[2] + " " + j_info[3]
			relationship_list.append([i, j_info[0]])

		count += 1 
		
	return relationship_list



def subtrees_function(root, subtrees = None):
	'''A function to check whether a node is one of various kinds,
	allowing us to check for duplications and hence make the largest
	possible subtrees without internal duplication events. Takes the 
	root node of the tree to be split into subtrees and returns a list
	of the subtrees. NB the subtree nodes in the list are still the 
	ones in the original tree, this doesn't make copies of them.'''
	
	if subtrees is None:
		subtrees = []

	children = [i for i in root.children]
	
	side1 = children[0]
	side2 = children[1]

	bipart1 = postorder3(side1)
	bipart2 = postorder3(side2)

	# Checking if *this* node corresponds to a duplication
	duplication = "No"

	for i in bipart1:
		if i in bipart2:
			duplication = "Yes"
	for i in bipart2:
		if i in bipart1:
			duplication = "Yes"

	# Checking if either side of the node contains another, downstream subduplication
	side1_dup = len(bipart1) - len(set(bipart1))
	side2_dup = len(bipart2) - len(set(bipart2))


	# If there are subduplications on both sides we need to go further into the tree
	if side1_dup != 0 and side2_dup != 0:
		subtrees_function(side1, subtrees)
		subtrees_function(side2, subtrees) 

	# If only one of the sides contains a subduplication, we make a subtree for the other side
	# and go further into the tree for this side
	elif side1_dup != 0:
		subtrees.append(side2)
		subtrees_function(side1, subtrees)
		
	elif side2_dup != 0:
		subtrees.append(side1)
		subtrees_function(side2, subtrees)
	
	# If there are no subduplications at all, great! 
	else:
		# If this node is a duplication node (it should be if we've 
		#got to this point), we make two subtrees
		if duplication == "Yes":
			subtrees.append(side1)
			subtrees.append(side2)
			
	return subtrees		
	

def tree_map(root, list):
	'''This replaces the labels of each node in a species tree with 
	the numbers of conflicts and concordances various gene trees show
	at that node. Takes the root of the tree to be mapped onto and a 
	list of bipartitions. WARNING: this function *replaces* the labels
	in the provided tree, it doesn't make a new tree with different ones
	new ones. Be careful when you call it.'''
	
	# Create a dictionary which records how many of each unique bipartition
	# there are in the list
	bipart_dict = {}

	for i in list:
		key = str(i[0][3:])
		if key in bipart_dict.keys():
			bipart_dict[key] += 1
		else:
			bipart_dict[key] = 1

	# Find each bipartition in the dictionary in the species tree using node_finder
	# and change the label to the number of times it was recorded in the list (i.e. 
	# the number of conflicts/concordances at that node in the gene tree(s))
	for i in bipart_dict.keys():	
		bipartition = i
		species_tree = root

		label = str(bipart_dict[i])
		node_finder(species_tree, bipartition, label)

def clear_labels(root):
	'''removes all the labels downstream of the root node specified, except
	taxon names'''

	for i in root.children:
		if i.istip == False:
			i.label = ''

		clear_labels(i)

def node_finder(root, bipartition, label):
	'''traverses tree and changes the label of the node with an identical bipartition 
	to the one specified'''
	
	for i in root.children:
		if i.istip == False:
			test_bipart = postorder3(i)
			test_bipart = str(test_bipart) # :/

			if test_bipart == bipartition:
				i.label = label
				
		node_finder(i, bipartition, label)


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python " + sys.argv[0] + " species_tree folder_of_homologs cutoff"
		sys.exit(0)
	
	# Initialise things
	cutoff = int(sys.argv[3])
	concordances = []
	conflicts = []
	
	# Open the first tree (species tree) file, build a tree and make bipartitions for it
	t1 = open(sys.argv[1],"r")
	for i in t1:
		n1 = i
	tree1, name_array1 = build(n1)
	tree1_biparts = postorder2(tree1)
	
	# Open and read the provided folder of homologs, go through each file in turn
	homologs_folder = sys.argv[2]
	file_list = os.listdir(sys.argv[2])
	
	for i in file_list:
		
		# Initialise
		relationship_list = []

		# Build tree 
		file_location = str(homologs_folder) + "/" + i
		t2 = open(file_location, "r")
		for i in t2:
			n2 = i
		tree2, name_array2 = build(n2)
		
		# Make subtrees and loop through computing biparts
		trees = subtrees_function(tree2)

		for i in trees:
			total_array = []
			i_biparts = postorder2(i, total_array, subtrees = True)
			print "\n\nnew tree: " + str(i_biparts)
			relationships = comp_biparts(tree1_biparts, i_biparts,name_array1, name_array2, sys.argv[2], cutoff)
			
			# Add each new bipart to the lists of concordances or conflicts
			for i in relationships:
				if i[1] == 'conflict':
					conflicts.append(i)
				elif i [1] == 'concordant':
					concordances.append(i)
		
	# Map concordances and conflicts back onto the tree using the lists
	print "\n"
	clear_labels(tree1)
	tree_map(tree1, concordances)
	concordance_tree = tree1.get_newick_repr(showbl = True)
	print concordance_tree + ";"

	clear_labels(tree1)
	tree_map(tree1, conflicts)
	conflict_tree = tree1.get_newick_repr(showbl = True)
	print conflict_tree + ";"
