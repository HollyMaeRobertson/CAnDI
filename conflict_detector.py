'''
This program is written to find gene tree conflict and concordance
between a species tree and a set of gene trees (homolog) trees, 
accounting for duplication events. When reverse is True it only takes
one homolog tree.

It can then return either (when 'normal') two species 
trees (in newick format) with node labels corresponding to the 
number of concordances or  conflicts at that node, 

or (when 'reverse') a single homolog tree with node labels 
showing conflict with the species tree, concordance with the species tree, 
or a duplication event.
'''

import sys, os
from node import Node
 
def build(instr):
	'''This takes in the tree in a newick format, as a string, then 
	puts it in a data structure that can be easily traversed'''

	root = None
	name_array = []
	index = 0
	nextchar = instr[index]
	begining = "Yep"
	keepgoing = True
	current_node = None
	counter = 0

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
			current_node.unique_id = str(counter)
			counter += 1
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
	

	return root, name_array


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

def change_tips_to_locus(root):
	'''postorder traversal of tree from specified root, changing
	all tips into unique values in the form species@locus'''
	
	for i in root.children:
		if i.istip:
			i.label = i.label + i.locus
		change_tips_to_locus(i)


def postorder2(root, total_list = None, subtrees = False):
	'''This traverses the entire tree downstream of the specified root
	and makes a list (total_list) of all the possible bipartitions, using postorder3.'''
 	
	# total_list should always start off empty
	if total_list is None:
		total_list = []

	# When subtrees is True, we want to make an extra bipartition that encompasses the whole subtree
	if subtrees == True:
		bipart = []
		other_side = [] # There is no other side as this is the whole tree
		
		bipart.append(str(root.unique_id))
		bipart.append(str(root.label))
		bipart.append(other_side)

		bipart = postorder3(root, bipart)

		total_list.append(bipart)
	
		subtrees = False

	# Making a bipart to add to total_list
	for i in root.children:
		if i.children:
			bipart = []

			bipart.append(str(i.unique_id)) 
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
	'''This function takes two biparts made by postorder2 and compares
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



def sort_conflicts(conflicts):
	'''This should sort a list of biparts that conflict 
	with the species tree into nodes and then from there 
	into separate "conflict topologies". Returns a
	dictionary which contains another dictionary for each
	node, which in turn contains the different conflicts.'''

	node_dict = {}

	for conflict in conflicts:
		bipart = conflict[2]
		
		bipart_node = conflict[0][0]

		if bipart_node not in node_dict.keys():
			node_dict[bipart_node] = [bipart]
		else:
			node_dict[bipart_node].append(bipart)

	complete_dict = {}


	for key in node_dict.keys():
		# List of all the biparts
		bipart_list = node_dict[key]
		
		# Make a 'holding' dictionary to put everything from this node into
		conflict_dict = {}

		# Add a first entry to the dictionary
		conflict_dict['conflict_1'] = [bipart_list[0]]

		# Now comes the hard bit
		for bipart in bipart_list[1:]:
			index = 1
			bp1 = bipart[3:]

			while True:

				name = "conflict_" + str(index)

				if name in conflict_dict.keys():
					bp2 = conflict_dict[name][0][3:]

				else:
					sys.stderr.write("Making list for node_" + key + " conflict_" + str(index) + "\r")
					conflict_dict[name] = [bipart]
					break

				rel = bipart_properties(bp1, bp2)
				
				if rel == 'concordant':
					conflict_dict[name].append(bipart)
					index = 1
					break
				
				else:
					# If not concordant, check next list entry
					index += 1

		complete_dict[key] = conflict_dict

	return complete_dict

def an_extremely_specific_function(a_list):
	'''This function takes a list and returns the length of 
	the second entry, as long as the second entry is also a list.
	It is mainly just for use in the below function conflict_stats, 
	and is hence  highly specific. '''
	if isinstance(a_list[1], list):
		return len(a_list[1])
	else:
		print "not a list!"

def conflict_stats(conflicts_dict, tree):
	'''This function should take a dictionary from sort_conflicts
	and calculate the most common conflict at each node, second-most
	common, etc.'''
	
	# First thing we're going to do is convert every entry in the 
	# dictionary within a dictionary to a list in form [key, [bp1, bp2, ...]]
	# Because *reasons*
	# (Maybe clean this up later Holly :/)
	stats_dict = {}

	for node in conflicts_dict.keys():
		stats_dict[node] = []
		for conflict in conflicts_dict[node].keys():
			bipart_list = conflicts_dict[node][conflict]
			new_list = [conflict, bipart_list]
			stats_dict[node].append(new_list)
		
	# Next we want to order all the conflicts within each node from most common to least common
	for node in stats_dict.keys():
		print ""
		node_on_tree = node_finder(tree, node)
		node_bipart = postorder3(node_on_tree)
		print "Node " + node + ": " + str(node_bipart) + "\n"
		stats_dict[node].sort(reverse=True, key=an_extremely_specific_function)

		# Then add up the total number of conflicts at that node? (where to store?)
		total = 0
		for conflict in stats_dict[node]:
			total += len(conflict[1])
		
		# Now let's calculate the percentages of conflict up to a certain threshold
		# Should threshold be % of total conflict or numbers of conflicts counted?
		counter = 0
		cumulative_percent = 0

		for conflict in stats_dict[node]:
			how_common = len(conflict[1])
			percent = float(how_common)/total * 100
			percent = round(percent, 2)
			
			if int(cumulative_percent) < 90 and counter < 5:
				print "At node " + str(node) + ", " + conflict[0] + " accounts for " + str(percent) + "% of all (" + str(total) + ") conflicts. The bipart is: " + str(conflict[1][0][3:])
			
			else:
				print "\nPercentage reported for node " + node + " (conflict_detector stops after 5 most common nodes currently): " + str(cumulative_percent) + "%. The remaining conflict is shared by " + str(len(stats_dict[node])-5) + " different conflicts.\n"
				break
			
			cumulative_percent += percent
			counter += 1



def make_string(bipart):
	'''This is to make it easier to understand the output when 
	printing things. It makes a neat string out of a bipart 
	made in postorder2 that shows both sides.'''
	
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
	
	NB: tree1 and tree2 are NOT equivalent. tree1 should be the tree you want to map your 
	relationships back onto.'''

	# Initialise things
	relationship_list = []
	test_bp1 = []
	test_bp2 = []
	rel = ''
	outf = open(log_name + ".log", "w")
	count = 0

	# Getting missing taxa so we can exclude them
	all_names = name_array1 + name_array2

	mis1, mis2 = unique_array(all_names, name_array1, name_array2)
	
	# Comparing all the biparts pairwise
	for i in tree1:
		various_relationships = []	
		lengths = []
		
		# Make the first bipart to test, removing missing taxa
		test_bp1 = list(set(i[3:]) - set(mis2)) 
		
		# When we remove missing taxa, we create a lot of duplicate bipartitions
		# We need to remove all but the most 'downstream' one
		# The most downstream one *should* be the one that was added to the list last...

		remaining_bps = [list(set(k[3:]) - set(mis2)) for k in tree1[count + 1:]]
	
		# So if there is another identical bipart further on in the list, we ignore this one
		if test_bp1 in remaining_bps:
			count += 1
			continue

		for j in tree2:
			
			# Make the second bipart to test
			test_bp2 = list(set(j[3:]) - set(mis1))

			# Mostly to make output legible
			not_bp1 = i[2]
			not_bp2 = j[2]
			bp1_string = make_string(i)			
			bp2_string = make_string(j)
			
			# What is the relationship between them?
			
			# We only want to record the relationship if the bootstrap > cutoff, OR if one or the other node doesn't have a bootstrap
			# Problematic as the label might not actually be a bootstrap even if it is an integer
			if str.isdigit(i[1]):
				i_cutoff = int(i[1])
			else:
				i_cutoff = 100

			if str.isdigit(j[1]):
				j_cutoff = int(j[1])
			else:
				j_cutoff = 100


			if i_cutoff < cutoff+1 or j_cutoff < cutoff+1:
				pass
			else:
				
				rel = bipart_properties(test_bp1, test_bp2)
				outf.write(str(rel) + ": " + bp1_string + i[1] + " " + bp2_string + j[1] + "\n")
				
				if rel == "conflict" or rel == "concordant":
					j_info = [str(rel), j]
					various_relationships.append(j_info)
					lengths.append(len(test_bp2))
					
					
			# Otherwise writes to the log file but doesn't store the relationship in the program
				else:
					outf.write(str(rel) + ": " + str(test_bp1) + i[1] +  " " + str(test_bp2) + j[1] + " "  + "Node not supported" + "\n")
		
		 
		# If there is only one thing in various_relationships, great! 
		# We put it and the node into the relationship_list
		if len(various_relationships) == 1:
			j_info = various_relationships[0]
			#print "\n" + str(count) + " " + j_info[0] + "\t" + bp1_string + " " + str(i[1]) + " " + str(i[0]) + "    " + j_info[1] + " " + j_info[2] + " " + j_info[3]
			relationship_list.append([i, j_info[0], j_info[1]])
			
		elif len(various_relationships) != 0:
			num = 100000000000000000
			counter = 0

			for k in lengths:
				if k < num:
					num = k
					j_info = various_relationships[counter]
				counter += 1
			
			#print "\n" + str(count) + " " + j_info[0] + "\t" + bp1_string + " " + str(i[1]) + " " + str(i[0]) + "    " + j_info[1] + " " + j_info[2] + " " + j_info[3]
			relationship_list.append([i, j_info[0], j_info[1]])

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
		if duplication == "Yes":
			subtrees.append(side1)
			subtrees.append(side2)
		elif duplication == "No":
			subtrees.append(root)
			
	return subtrees		


def label_duplications(root, recursive = True):
	'''Labels duplication nodes with "D"'''
	if root.istip == False:
		children = [i for i in root.children]
	
		side1 = children[0]
		side2 = children[1]
		
		if side1.istip:
			bipart1 = [side1.label]
		else:
			bipart1 = postorder3(side1)
		
		if side2.istip:
			bipart2 = [side2.label]
		else:
			bipart2 = postorder3(side2)
		
		duplication = "No"
		
		for i in bipart1:
			if i in bipart2:
				duplication = "Yes"
		for i in bipart2:
			if i in bipart1:
				duplication = "Yes"

		if duplication == "Yes":	
			root.label = 'D'
	if recursive == True:
		for node in root.children:
			label_duplications(node)
		

def tree_map(root, bipart_list):
	'''This replaces the labels of each node in a species tree with 
	the numbers of conflicts and concordances various gene trees show
	at that node. Takes the root of the tree to be mapped onto and a 
	list of bipartitions. 
	
	WARNING: this function *replaces* the labels
	in the provided tree, it doesn't make a new tree with different ones
	new ones. Be careful when you call it.'''
	
	# Create a dictionary which records how many of each unique bipartition
	# there are in the list
	bipart_dict = {}

	for i in bipart_list:
		key = str(i[0][0])
		if key in bipart_dict.keys():
			bipart_dict[key] += 1
		else:
			bipart_dict[key] = 1
	
	# Find each bipartition in the dictionary in the species tree using node_finder
	# and change the label to the number of times it was recorded in the list (i.e. 
	# the number of conflicts/concordances at that node in the gene tree(s))
	for key in bipart_dict.keys():
		label = str(bipart_dict[key])
		node = node_finder(root, key)
		node.label = label


def tree_map2(root, rel_list, label):
	'''this replaces the labels of each node in a species tree with a given label,
	provided they are in the list of biparts passed to the function (rel_list)'''

	bipart_dict = {}
	if root.parent:
		root = root.parent

	for i in rel_list:
		key = str(i[0][0])
		bipart_dict[key] = label
	

	for key in bipart_dict.keys():
		node = node_finder(root, key)
		node.label = label

def clear_labels(root):
	'''removes all the labels downstream of the root node specified, except
	taxon names'''
	if root.istip == False:
		root.label = ''

	for i in root.children:
		if i.istip == False:
			i.label = ''

		clear_labels(i)

def change_tips_to_species(root):
	'''replaces names in the form species@locus with only the species name'''

	for i in root.children:
		if i.istip == True:
			name = ''
			for j in range(len(i.label)):
				if i.label[j] != '@':
					name += i.label[j]
				else:
					i.label = name
					break
		
		change_tips_to_species(i)	

def node_finder(root, node_id):
	'''Calls node_finder_recursive to find the node, then returns it.
	Bit of a hack but ok'''
	
	node = node_finder_recursive(root, node_id, node = None)
	if len(node) == 1:
		node = node[0]
	else:
		print "problem :("

	return node


def node_finder_recursive(root, node_id, node = None, first_time = True):
	'''Traverses tree and returns the node with the corresponding id'''

	if node is None:
		node = []
	
	if first_time == True:
		if root.unique_id == node_id:
			node.append(root)

	for i in root.children:
		if i.istip == False:
			current_id = i.unique_id
			
			if current_id == node_id:
				node.append(i)
								
		node_finder_recursive(i, node_id, node, first_time = False)

	return node
		

def compare_trees(tree1_biparts, name_array1, tree2, mode, cutoff):
	'''This function compares a subtree to tree1, which has already
	been split into biparts.'''

	conflicts = []
	concordances = []
	
	# name_array2 must also include the species from the next species up (even though 
	# nothing should really be compared against them except for the subtree root)
	keepgoing = True
	current_node = tree2
	new_names = []
	
	while True:
		parent = current_node.parent
		if parent != None:	
			label_duplications(parent, recursive = False)
			if parent.label == 'D':
				current_node = parent
			else:
				for i in parent.children:
					if i != current_node:
						new_names = postorder3(i)
				break
		else:
			break
		
	name_array2 = postorder3(tree2)
	name_array2.extend(new_names)

	tree2_biparts = postorder2(tree2, subtrees = True)
		
	if mode == "normal":
		relationships = comp_biparts(tree1_biparts, tree2_biparts, name_array1, name_array2, sys.argv[2], cutoff)
	elif mode == "reverse":
		relationships = comp_biparts(tree2_biparts, tree1_biparts, name_array2, name_array1, sys.argv[2], cutoff)
	
	for relationship in relationships:
		if relationship[1] == 'conflict':
			conflicts.append(relationship)
		elif relationship[1] == 'concordant':
			concordances.append(relationship)

	return conflicts, concordances

def identify_tricky_nodes(root, subtree_list, tricky_nodes = None, is_root = True):
	'''This takes the root node of a tree, a list of all the subtree nodes, 
	and returns a list of all the nodes that are not either a) a duplication
	node; b) a subtree root or c) within a subtree '''
	
	if tricky_nodes is None:
		tricky_nodes = []
	
	if is_root == True:
		label_duplications(root, recursive = False)

		if root.label != 'D' and root not in subtree_list and root.istip == False:
			tricky_nodes.append(root)
		elif root in subtree_list:
			return tricky_nodes
	
	for child in root.children:
		label_duplications(child, recursive = False)
		
		if child.label != 'D' and child not in subtree_list and child.istip == False:
			tricky_nodes.append(child)
		elif child in subtree_list or child.istip == True:
			continue
		identify_tricky_nodes(child, subtree_list, tricky_nodes, is_root = False)

	return tricky_nodes

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python " + sys.argv[0] + " species_tree folder_of_homologs cutoff mode(normal/reverse/bipart_to_search in form X,Y,Z)"
		sys.exit(0)
	
	# Initialise things
	cutoff = int(sys.argv[3])
	mode = sys.argv[4]
	
	# Open the first tree (species tree) file, build a tree and make bipartitions for it
	t1 = open(sys.argv[1],"r")
	for i in t1:
		n1 = i
	tree1, name_array1 = build(n1)
	tree1_biparts = postorder2(tree1)
	all_taxa = []
	all_taxa.append(tree1.unique_id)
	all_taxa.append(tree1.label)
	all_taxa.append([])
	all_taxa = postorder3(tree1, all_taxa)	
	tree1_biparts.append(all_taxa)
	
	# When we're reading a folder of homologous gene trees
	if mode == "normal":
		# Initialise
		total_conflicts = []
		total_concordances = []

		# Open and read the provided folder of homologs, go through each file in turn
		homologs_folder = sys.argv[2]
		file_list = os.listdir(sys.argv[2])
	

		for file in file_list:
			sys.stderr.write("processing " + str(file) + '\r')

			# Build tree 
			file_location = str(homologs_folder) + "/" + str(file)
			t2 = open(file_location, "r")
			for i in t2:
				n2 = i
			tree2, name_array2 = build(n2)
			

			# Make subtrees
			trees = subtrees_function(tree2)

			# Loop through finding conflicts and concordances

			for tree in trees:
				conflicts, concordances = compare_trees(tree1_biparts, name_array1, tree, mode, cutoff)
				total_conflicts.extend(conflicts)
				total_concordances.extend(concordances)
				
			# Find the tricky nodes
			tricky_nodes = identify_tricky_nodes(tree2, trees)
			
		
			
			for node in tricky_nodes:
				node_bipart = []
				node_bipart_list = []
				node_bipart.append(node.unique_id)
				node_bipart.append(node.label)
				node_bipart.append([])
				node_bipart = postorder3(node, node_bipart)
				node_bipart_list.append(node_bipart)
				name_array2 = postorder3(node)
				new_names = []

				current_node = node

				while True:
					parent = current_node.parent
					if parent == None:
						break
					
					label_duplications(parent, recursive = False)
					if parent.label == 'D':
						current_node = parent
					else:
						for i in parent.children:
							if i != current_node:
								new_names = postorder3(i)
						break
				
				name_array2.extend(new_names)

				rel_list = comp_biparts(tree1_biparts, node_bipart_list, name_array1, name_array2, sys.argv[2], cutoff)
				
				for rel in rel_list:
						
					if rel[1] == 'conflict':
						total_conflicts.append(rel)
					elif rel[1] == 'concordant':
						total_concordances.append(rel)
				
		total_conflicts_1 = total_conflicts[:]
		sorted_conflicts = sort_conflicts(total_conflicts_1)
		conflict_stats(sorted_conflicts, tree1)

		# Map concordances and conflicts back onto the tree using the lists
		print "\n"
		clear_labels(tree1)
		tree_map(tree1, total_concordances)
		concordance_tree = tree1.get_newick_repr(showbl = True)
		print "Concordance tree: "
		print concordance_tree + ";"

		clear_labels(tree1)
		tree_map(tree1, total_conflicts)
		conflict_tree = tree1.get_newick_repr(showbl = True)
		print "Conflict tree:"
		print conflict_tree + ";"

	# When we're reading one homologous gene tree and mapping conflict/concordance with the species tree back onto that
	elif mode == 'reverse':
		
		# Open the homologous gene tree file and build a tree
		t2 = open(sys.argv[2], "r")
		for i in t2:
			n2 = i
		tree2, name_array2 = build(n2)

		# Make subtrees
		trees = subtrees_function(tree2)
		
		if len(trees) != 0:

			# Finding conflicts and concordances and mapping onto the species tree
			clear_labels(tree2)
			
			for tree in trees:
				conflicts, concordances = compare_trees(tree1_biparts, name_array1, tree, mode, cutoff)
				tree_map2(tree, conflicts, 'X')
			
			
			for conflict in conflicts:
				print conflict

			for tree in trees:
				conflicts, concordances = compare_trees(tree1_biparts, name_array1, tree, mode, cutoff)
				tree_map2(tree, concordances, '*')
			
			# Find the tricky nodes
			tricky_nodes = identify_tricky_nodes(tree2, trees)
			conflicts = []
			concordances = []

			for node in tricky_nodes:
				node_bipart = []
				node_bipart_list = []
				node_bipart.append(node.unique_id)
				node_bipart.append(node.label)
				node_bipart.append([])
				node_bipart = postorder3(node, node_bipart)
			
				node_bipart_list.append(node_bipart)
				name_array2 = postorder3(node)
				new_names = []	
				current_node = node

				while True:
					parent = current_node.parent
					if parent != None:
						label_duplications(parent, recursive = False)
						if parent.label == 'D':
							current_node = parent
						elif type(parent) is None:
							break
						else:
							for i in parent.children:
								if i != current_node:
									new_names = postorder3(i)
							break
					else:
						break
					
				name_array2.extend(new_names)

				rel_list = comp_biparts(node_bipart_list, tree1_biparts, name_array2, name_array1, sys.argv[2], cutoff)

				for rel in rel_list:

					if rel[1] == 'conflict':
						conflicts.append(rel)
					elif rel[1] == 'concordant':
						concordances.append(rel)
			
			#map these back on to whole tree
			print conflicts

			tree_map2(tree2, conflicts, 'X')
			tree_map2(tree2, concordances, '*')

			label_duplications(tree2)

			change_tips_to_locus(tree2)
			new_tree = tree2.get_newick_repr(showbl = True)
			print new_tree + ";"
		
		else:
			conflicts, concordances = compare_trees(tree1_biparts, name_array1, tree2, mode, cutoff)

			#map these back on to whole tree
			print conflicts

			tree_map2(tree2, conflicts, 'X')
			tree_map2(tree2, concordances, '*')

			label_duplications(tree2)

			change_tips_to_locus(tree2)
			new_tree = tree2.get_newick_repr(showbl = True)
			print new_tree + ";"


	else:
		
		# Initialise
		total_conflicts = []
		total_concordances = []

		# Put the query bipart in a form comp_biparts likes
		query_entered = mode.split(",")
		query = [['', '', []]]
		for taxon in query_entered:
			query[0].append(taxon)

		print query

		query_name_array = query_entered
		print query_name_array

		# Open and read the provided folder of homologs
		homologs_folder = sys.argv[2]
		file_list = os.listdir(sys.argv[2])
	
		concordant_files = []

		for file in file_list:
			sys.stderr.write("processing " + str(file) + '\r')

			# Build tree 
			file_location = str(homologs_folder) + "/" + str(file)
			t2 = open(file_location, "r")
			for i in t2:
				n2 = i
			tree2, name_array2 = build(n2)

			# Make subtrees
			trees = subtrees_function(tree2)

			# Look for concordances in the subtrees
			for tree in trees:
				conflicts, concordances = compare_trees(query, name_array1 , tree, 'normal', cutoff)
			
			# Look for concordances in the tricky nodes 
			# NB This is copy-pasted code >:-(
			# Consider putting it in a function when you clean it up, Holly
				
			tricky_nodes = identify_tricky_nodes(tree2, trees)
			
			for node in tricky_nodes:
				node_bipart = []
				node_bipart_list = []
				node_bipart.append(node.unique_id)
				node_bipart.append(node.label)
				node_bipart.append([])
				node_bipart = postorder3(node, node_bipart)
				node_bipart_list.append(node_bipart)
				name_array2 = postorder3(node)
				new_names = []

				current_node = node

				while True:
					parent = current_node.parent
					if parent == None:
						break
					
					label_duplications(parent, recursive = False)
					if parent.label == 'D':
						current_node = parent
					else:
						for i in parent.children:
							if i != current_node:
								new_names = postorder3(i)
						break
				
				name_array2.extend(new_names)

				rel_list = comp_biparts(query, node_bipart_list, name_array1, name_array2, sys.argv[2], cutoff)
				
				for rel in rel_list:			
					if rel[1] == 'concordant':
						concordances.append(rel)
			
			#Check concordance list for concordances
			if len(concordances) != 0:
				concordant_files.append(file)
			

		
		#for file in concordant_files:
			#print str(file)
		
		print "There are " + str(len(concordant_files)) + " files containing the specified bipart, out of " + str(len(file_list)) + " files in total."

