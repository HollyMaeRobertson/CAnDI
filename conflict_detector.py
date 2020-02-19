'''
Program to get tree length, concordance or whatever
'''
import sys
from node import Node
 
#This takes in the newick and the
#seq data then puts them in a data
#structure that can be preorder or
#postorder traversed pretty easily
def build(instr):
	#print "Entered build"
	root = None
	name_array =[]
	index = 0
	nextchar = instr[index]
	begining = "Yep"
	keepgoing = True
	current_node = None
	#keeps going until the value becomes false
	while keepgoing == True:
		#This situation will only happen at the very beginning but
		#when it hits this it will create a root and change begining
		#to no
		if nextchar == "(" and begining == "Yep":
				
			root = Node()
			current_node = root
			begining = "No"
		#This happens anytime their is an open bracket thats not the
		#beginning
		elif nextchar == "(" and begining == "No":
		
			newnode = Node()
			current_node.add_child(newnode)
			current_node = newnode
		#This indicates that you are in a clade and tells the 
		#program to move back one to grab the sister to the clade
		elif nextchar == ',':
		
			current_node = current_node.parent
		#This says you are closing a clade and therefore it moves
		#back to where the parent node is which allows the name
		#to be added to the parent node
		elif nextchar == ")":
			#print "Closing Clade"
			current_node = current_node.parent
			index += 1
			nextchar = instr[index]
			#this puts bootstrap labels on nodes
			while True:
			
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			current_node.label = name
			index -= 1
		#This indicates everything is done so keepgoing becomes false
		elif nextchar == ';':
		
			keepgoing = False
			break
		#This indicates you have branch lengths so it grabs the branch
		#lengths turns them into floats and puts them in the current node
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
		#This is for if anywhitespace exists
		elif nextchar == ' ':
		
			pass #means the whitespace doesn't do anything weird
            	
		#this is for if there are locus labels
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
				
		
		#This is for when any taxa name is hit, it will concatenate
		#the taxa names together and add the name
		else: # this is an external named node
		
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
		if index < len(instr) - 1:
			index += 1
		nextchar = instr[index]
		name = ""
		branch = ""
	return root,name_array


def postorder3(root,bipart):
	for i in root.children:
		if i.istip:
			bipart.append(i.label)
	
		postorder3(i,bipart) #---> recursion

def postorder2(root,total_array, subtrees):

	cutoff = 0
 	
	if subtrees == True:
		#all tips go into one bipartition
		bipart = []
		bipart.append(str(root.length))
		bipart.append(str(root.label))
	
		other_side = []

		bipart.append(other_side)

		postorder3(root, bipart)
		
		total_array.append(bipart)
	
	subtrees = False

	for i in root.children:
	#part has children so grab'em
		if i.children:
			#checks internal labels
			#if i.label 
			bipart = []
			bipart.append(str(i.length)) #gives length of node
			bipart.append(str(i.label)) #e.g. bootstrap values
			
			#what is on the other side?
			other_side = []

			'''we want the children etc of i.parent that 
			aren't the children of i - to go in not_bipart'''
			one_up = i.parent
			for j in one_up.children:
				if j == i:
					pass
				else:
					if j.istip:
						other_side.append(j.label)
					postorder3(j, other_side)
				
			#we want bipart to include the 'not_bipart' information	
			bipart.append(other_side)
		
			#the actual bipart
			postorder3(i,bipart)

			total_array.append(bipart)


		postorder2(i,total_array, subtrees = False)  #---> recursion
	
		
	return total_array



def unique_array(array,tree1,tree2):
	all_species = set()
	for x in array:
		all_species.add(x)
	mis1 = list(set(all_species) - set(tree1)) #the species that are in 2 but not 1
	mis2 = list(set(all_species) - set(tree2)) #the species that are in 1 but not 2
	return mis1,mis2

'''
How do the biparts relate
'''
def bipart_properties(bp1,bp2):
	
	keepgoing = "false"
	#difference in the biparts
	diff1 = list(set(bp1) - set(bp2))
	diff2 = list(set(bp2) - set(bp1))

	#print "bipart 1: " + str(bp1)
	#print "bipart 2: " + str(bp2)
	#print diff1
	#print diff2
	#check for no overlp
	if len(bp1) == 1 or len(bp2) == 1:
		#print "uninformative"
		return "uninformative"
	elif len(diff1) == len(bp1):
		#print "no comp"
		return "no_comp"
	#check if nested
	elif len(bp1) == len(bp2) and len(diff1) == 0:
		#print "identical"
		return "concordant"
	elif len(diff1) == 0:
		#print "nested in 2"
		return "1 nested in 2"
	elif len(diff2) == 0:
		#print "nested in 1"
		return "2 nested in 1"
	else:
		#print "conflict"
		return "conflict"

	

def make_string(bipart, not_bipart):
	'''purely for aestheics'''

	name = ''
	for i in bipart:
		name += i + " "
	name += '| '
	for i in not_bipart:
		name += i + " "

	return str(name)

'''
Compare the bipartitions
'''
def comp_biparts(tree1,tree2,name_array1,name_array2,log_name,cutoff):
	
	all_names = []
	test_bp1 = []
	test_bp2 = []
	rel = 0
	all_species = set()
	outf = open(log_name + ".log", "w")
	
	'''
	get names to know what can and can't be evaluated
	'''

	all_names = name_array1 + name_array2
	# get the names, for missing taxa check etc.
	mis1,mis2 = unique_array(all_names,name_array1,name_array2)
	count = 0
	for i in tree1:
		
		#list for the things to go into
		various_relationships = []	
		lengths = []

		for j in tree2:
			'''i is bipart from tree1 (species tree), check for concordance
			unable to speak to, and conflict for everything in tree2
			remove stuff that's not going to be found first'''		
		
			test_bp1 = list(set(i[3:]) - set(mis2)) 
			test_bp2 = list(set(j[3:]) - set(mis1))
			
					
			#the other side
			not_bp1 = i[2]
			not_bp2 = j[2]

			#make a string to print (just to make things look nicer)
			bp1_string = make_string(test_bp1, not_bp1)
			bp2_string = make_string(test_bp2, not_bp2)
			
			rel = bipart_properties(test_bp1, test_bp2)

			#compare them if bootstrap > cutoff
			#if there is no boostrap, bootstrap defaults to 
			if type(i[1]) == str or type(j[1]) == str or str(i[1]) == "" or str(j[1]) == ""  or int(i[1]) > cutoff and int(j[1]) > cutoff: 
				outf.write(str(rel) + ": " + str(test_bp1) + i[1] + " " + str(test_bp2) + j[1] + "\n")
				
				if rel == "conflict" or rel == "concordant":
					j_info = [str(rel), bp2_string, str(j[1]), str(j[0])]
					various_relationships.append(j_info)
					lengths.append(len(test_bp2))
					
			#otherwise write out that the node is unsupported
			else:
				outf.write(str(rel) + ": " + str(test_bp1) + i[1] +  " " + str(test_bp2) + j[1] + " "  + "Node not supported" + "\n")
		
		
		if len(various_relationships) == 1:
			j_info = various_relationships[0]
			print "\n" + str(count) + " " + j_info[0] + "\t" + bp1_string + " " + str(i[1]) + " " + str(i[0]) + "    " + j_info[1] + " " + j_info[2] + " " + j_info[3]
		elif len(various_relationships) != 0:

			num = 9900000000000000000009
			counter = 0

			for k in lengths:
				if k < num:
					num = k
					j_info = various_relationships[counter]
				counter += 1
			
			print "\n" + str(count) + " " + j_info[0] + "\t" + bp1_string + " " + str(i[1]) + " " + str(i[0]) + "    " + j_info[1] + " " + j_info[2] + " " + j_info[3]
			
		count += 1 
		
		


'''OK: for each node, we want to look at everything either side of the bipart'''

def quick_biparts(root, node_bipart):

	for i in root.children:
		if i.istip:	
			node_bipart.append(i.label)
		
		quick_biparts(i, node_bipart)
			
	return node_bipart


def what_kind_of_node(root, subtrees):

	children = [i for i in root.children]
	
	side1 = children[0]
	side2 = children[1]

	bipart1 = quick_biparts(side1, [])
	bipart2 = quick_biparts(side2, [])

	#is there a duplication?
	duplication = 0 

	for i in bipart1:
		if i in bipart2:
			duplication += 1
	for i in bipart2:
		if i in bipart1:
			duplication += 1

	
	#there was a duplication
	#is there a subduplication on either side
	side1_dup = len(bipart1) - len(set(bipart1))
	side2_dup = len(bipart2) - len(set(bipart2))

	#case that both subduplication
	if side1_dup != 0 and side2_dup != 0:
		what_kind_of_node(side1, subtrees)
		what_kind_of_node(side2, subtrees) #recursion!

	#cases that one subduplication
	elif side1_dup != 0:
		subtrees.append(side2)
		what_kind_of_node(side1, subtrees)
		
	elif side2_dup != 0:
		subtrees.append(side1)
		what_kind_of_node(side2, subtrees)
	
	#case that no subduplications
	else:
		#if this node is a duplication node
		if duplication != 0:
			subtrees.append(side1)
			subtrees.append(side2)
			
	return subtrees		
	

def subtrees_function(root):

	#initialise a list of 'root' nodes of subtrees
	subtrees = []

	whatever = what_kind_of_node(root, subtrees)

	return whatever





if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python " + sys.argv[0] + " treefile1 treefile2 cutoff"
		sys.exit(0)
	
	tree1_biparts = []
	tree2_biparts = []
	total_array = []
	name_array1 = []
	name_array2 = []
	cutoff = int(sys.argv[3])

	t1 = open(sys.argv[1],"r")
	for i in t1:
		n1 = i
	tree1,name_array1 = build(n1)
	#at this point we want to split the tree into subtrees?
	
	tree1_biparts = postorder2(tree1, total_array, subtrees = False)
	
	t2 = open(sys.argv[2],"r")
	for i in t2:
		n2 = i
	tree2,name_array2 = build(n2)
	total_array = []
	tree2_biparts = postorder2(tree2, total_array, subtrees = False)
	
#	comp_biparts(tree1_biparts,tree2_biparts,name_array1,name_array2, sys.argv[2], cutoff)
	
	trees = subtrees_function(tree2)

	#
	for i in trees:
		total_array = []
		i_biparts = postorder2(i, total_array, subtrees = True)
		print "\n\nnew tree: " + str(i_biparts)
		comp_biparts(tree1_biparts, i_biparts,name_array1, name_array2, sys.argv[2], cutoff)

