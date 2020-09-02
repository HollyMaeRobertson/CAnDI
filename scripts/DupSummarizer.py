
"""
this takes a newick string as instr
and reads the string and makes the 
nodes and returns the root node
"""
#This code will let you walk through a tree
#and act on nodes that you want to
import sys
from node import Node
 
#This takes in the newick and the
#seq data then puts them in a data
#structure that can be preorder or
#postorder traversed pretty easily
def build(instr):
	#print "Entered build"
	root = None
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
		
			index += 1
			nextchar = instr[index]
		#This is for when any taxa name is hit, it will concatenate
		#the taxa names together and add the name
		else: # this is an external named node
		
			newnode = Node()
			current_node.add_child(newnode)
			current_node = newnode
			current_node.istip = True
			while True:
				if nextchar == ',' or nextchar == ')' or nextchar == ':' \
					or nextchar == ';' or nextchar == '[':
					break
				name += nextchar
				index += 1
				nextchar = instr[index]
			current_node.label = name
			index -= 1
		if index < len(instr) - 1:
			index += 1
		nextchar = instr[index]
		name = ""
		branch = ""
	return root


'''
Gets all the nodes (long story... Im an idiot) in a tree structure
'''
def postorder_name_getter(side,root,array):
	
	if len(array) == 0:
		array.append(root.label)
	for i in root.children:
		if i.istip == False:

			array.append(i.label)
        
		postorder_name_getter(side,i,array)

'''
Does the traversal from root to tip
'''
def preorder(root,HASH,tip_HASH,len_HASH,Parent_HASH,Species_HASH,counter,count_array):
	
	
	for child in root.children:

		if child.label == "D":
			
			aa = []
			tip_aa = []
			tlen_aa = []
			
			#count array provides a hacky way to not lose info during the recursive
			#process
			if len(count_array) != 0:
				new_count = int(count_array[-1]) + 1
				count_array.append(new_count)
			else:
				new_count = 1
				count_array.append(new_count)
			

			Parent_HASH[count_array[-1]] = child.length
			
			array = []
			tip_array = []
			tlen_array = []
			
			postorder_name_getter("left",child.children[0],array)
			aa.append(array)
			

			if child.children[0].istip:
				tlen_array.append(child.children[0].length)
				tip_array.append(child.children[0].label)
			else:
				tree_length_get(child.children[0],tlen_array)
				tip_name_get(child.children[0],tip_array)
			tlen_aa.append(tlen_array)
			tip_aa.append(tip_array)
			
			array = []
			tip_array = []
			tlen_array = []
			
			postorder_name_getter("right",child.children[1],array)
			aa.append(array)
			
			
			if child.children[1].istip:
				tlen_array.append(child.children[1].length)
				tip_array.append(child.children[1].label)
			else:
				tree_length_get(child.children[1],tlen_array)
				tip_name_get(child.children[1],tip_array)
			tlen_aa.append(tlen_array)
			tip_aa.append(tip_array)
			
			species_bipart = []
			species_bipart = check_corresponding_species_bipart(tip_aa,bipart_array)

			
			Species_HASH[count_array[-1]] = ",".join(species_bipart)
			
			HASH[count_array[-1]] = aa
			tip_HASH[count_array[-1]] = tip_aa
			len_HASH[count_array[-1]] = tlen_aa

		preorder(child,HASH,tip_HASH,len_HASH,Parent_HASH,Species_HASH,counter,count_array)

'''
Gets tip names for a subtree
'''
def tip_name_get(root,tip_array):
	
	for child in root.children:
		if child.istip:
			tip_array.append(child.label)
	
		tip_name_get(child,tip_array)

'''
Get the branches
'''
def tree_length_get(root,tlen_array):
	
	if root.length:
		tlen_array.append(root.length)
	for child in root.children:

		tree_length_get(child,tlen_array)		

'''
Get the length of edge that proceeds each thing (duplication,conflict,concordance)
'''
def preorder2(root,Association_HASH):

	for child in root.children:
		
		if child.istip == False:
			if child.label:
				Association_HASH[child.label].append(child.length)
				
		
		preorder2(child,Association_HASH)

'''
Takes subtree and gets tips
'''
def get_tips(root,array):

	for child in root.children:
		if child.istip:
			array.append(child.label)
		get_tips(child,array)

'''
Get the biparts associated with the species tree
'''
def get_biparts(root,bipart_array):
	
	for child in root.children:
		if child.istip == False:
			array = []
			get_tips(child,array)
			bipart_array.append(array)
		get_biparts(child,bipart_array)

'''
Get all taxa, semi taken from comparisons
'''
def all_unique(tree1, tree2):

	
	# Make a set of all the species in both trees, so we only have unique 
	# species.
	array = tree1 + tree2
	all_species = set()
	for x in array:
		
		all_species.add(x.split("@")[0])
	
	set_array = []
	for x in all_species:
		set_array.append(x)
	return set_array
			
'''
Check the bipart in species array associated with duplication
'''
def check_corresponding_species_bipart(dup_array,species_array):


	unique_array = []
	sort_array = []		
	unique_array = all_unique(dup_array[0],dup_array[1])
	
	#Need all equal taxa and get the smallest where they have all equal taxa
	for x in species_array:
		flag = 0
		if(set(unique_array).issubset(set(x))): 
   			 flag = 1
	
		if flag == 1:
			sort_array.append(x)
	
	sort_array.sort(key=len)
	if 1 < len(sort_array) and len(sort_array[0]) == len(sort_array[1]):
		print "Two locations for Dup?"
		sys.exit(0)
	
	return sort_array[0]



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("python "+sys.argv[0]+" DuplicationTree SpeciesTree")
        sys.exit(0)

name = sys.argv[1]
name2 = open(name, "r");

sp_t = sys.argv[2]
sp_t2 = open(sp_t,"r")

outw = open(name+".dupinfo.tsv","w")
outw2 = open(name+".treesum.tsv","w")

#Turn species tree into biparts
for i in sp_t2:
	n = build(i)
	bipart_array = []
	all_array_for_out = []
	get_tips(n,all_array_for_out)
	get_biparts(n,bipart_array)
	bipart_array.append(all_array_for_out)

for i in name2:
	n = build(i)
	#Contains node info
	HASH = {}
	#contains tip info
	tip_HASH = {}
	#edge length
	len_HASH = {}
	#Parent of Dup
	Parent_HASH = {}
	#Store part of species tree duplication is associated with
	Species_HASH = {}
	
	
	#Total Numbers
	Association_HASH = {}
	Association_HASH["D"] = []
	Association_HASH["*"] = []
	Association_HASH["X"] = []
	Association_HASH["U"] = []
	counter = 0
	count_array = []

	'''
	Check the very first node
	'''
	if n.label == "D":
		
		aa = []
		tip_aa = []
		counter += 1
		count_array.append(counter)
		array = []
		tip_array = []
		tree_length_array = []
		tlen_aa = []
		tlen_array = []
		Association_HASH["D"].append("NA")
		
		postorder_name_getter("left",n.children[0],array)
		aa.append(array)
		

		
		
		
		if n.children[0].istip:
			tip_array.append(n.children[0].label)
			tlen_array.append(n.children[0].length)
		else:
			tree_length_get(n.children[0],tlen_array)
			tip_name_get(n.children[0],tip_array)
		tlen_aa.append(tlen_array)
		tip_aa.append(tip_array)
		
		array = []
		tip_array = []
		tlen_array = []
		postorder_name_getter("right",n.children[1],array)
		aa.append(array)
		
		if n.children[1].istip:
			tip_array.append(n.children[1].label)
			tlen_array.append(n.children[1].length)
		else:
			tip_name_get(n.children[1],tip_array)
			tree_length_get(n.children[1],tlen_array)
		
		tip_aa.append(tip_array)
		tlen_aa.append(tlen_array)
		
		species_bipart = []
		species_bipart = check_corresponding_species_bipart(tip_aa,bipart_array)

		Species_HASH[counter] = ",".join(species_bipart)
		len_HASH[counter] = tlen_aa
		tip_HASH[counter] = tip_aa
		HASH[counter] = aa
	elif n.label == "*":
		Association_HASH["*"].append("NA")
	elif n.label == "X":
		Association_HASH["X"].append("NA")
	elif n.label == "U":
		Association_HASH["U"].append("NA")
		

	preorder(n,HASH,tip_HASH,len_HASH,Parent_HASH,Species_HASH,counter,count_array)

	#Summarize left and right
	outw.write("Dup#\tLengthOfEdge\tLeft_Duplications\tRight_Duplication\tDifferenceInDuplications\tLeft_Conflict\tRight_Conflict\tDifferenceInConflict\tLeft_Concordance\tRight_Concordance\tDifferenceInConcordance\tLeft_Uninformative\tRight_Uninformative\tDifferenceInUninformative\tLeftTips\tRightTips\tDifferenceInTips\tLeftAverageEdgeLength\tRightAverageEdgeLength\tDiffEdgeLength\tOrthologBipart\n")
	for x in HASH:

		Left_Duplication = 0
		Left_Conflict = 0
		Left_Concordance = 0
		Left_Uninformative = 0
		#left side
		for i in HASH[x][0]:
			if i == "D":
				Left_Duplication += 1
			elif i == "*":
				Left_Concordance += 1
			elif i == "X":
				Left_Conflict += 1
			elif i == "U":
				Left_Uninformative += 1
				
			
			Right_Duplication = 0
			Right_Conflict = 0
			Right_Concordance = 0
			Right_Uninformative = 0
			for i in HASH[x][1]:
				if i == "D":
					Right_Duplication += 1
				elif i == "*":
					Right_Concordance += 1
				elif i == "X":
					Right_Conflict += 1
				elif i == "U":
					Right_Uninformative += 1
				else:
					print "Somethings Wrong with node labels"
			 


		
		Left_Tips = len(tip_HASH[x][0])
		Right_Tips = len(tip_HASH[x][1])
		
		Left_Sum = (sum(len_HASH[x][0]) / len(len_HASH[x][0]))
		Right_Sum = (sum(len_HASH[x][1]) / len(len_HASH[x][1]))
		
		DupDiff = abs(Left_Duplication - Right_Duplication)
		ConfDiff = abs(Left_Conflict - Right_Conflict)
		ConcDiff = abs(Left_Concordance - Right_Concordance)
		UnDiff = abs(Left_Uninformative - Right_Uninformative)
		TipDiff = abs(Left_Tips - Right_Tips)
		SumDiff = abs(Left_Sum - Right_Sum)
		
		if x in Parent_HASH:
			outw.write("Duplication" + str(x) + "\t" + str(Parent_HASH[x]) + "\t" + str(Left_Duplication) + "\t" + str(Right_Duplication) + "\t" + str(DupDiff) \
			+ "\t" + str(Left_Conflict) + "\t" + str(Right_Conflict) + "\t" + str(ConfDiff) + "\t" \
			+ str(Left_Concordance) + "\t" + str(Right_Concordance) + "\t" + str(ConcDiff) + "\t" \
			+ str(Left_Uninformative) + "\t" + str(Right_Uninformative) + "\t" + str(UnDiff) + "\t" \
			+ str(Left_Tips) + "\t" + str(Right_Tips) + "\t" + str(TipDiff) + "\t" + str(Left_Sum) + "\t" \
			+ str(Right_Sum) + "\t" + str(SumDiff) + "\t" + Species_HASH[x] + "\n")

		else:
			outw.write("Duplication" + str(x) + "\t" + "NA" + "\t" + str(Left_Duplication) + "\t" + str(Right_Duplication) + "\t" + str(DupDiff) \
			+ "\t" + str(Left_Conflict) + "\t" + str(Right_Conflict) + "\t" + str(ConfDiff) + "\t" \
			+ str(Left_Concordance) + "\t" + str(Right_Concordance) + "\t" + str(ConcDiff) + "\t" \
			+ str(Left_Uninformative) + "\t" + str(Right_Uninformative) + "\t" + str(UnDiff) + "\t" \
			+ str(Left_Tips) + "\t" + str(Right_Tips) + "\t" + str(TipDiff) + "\t" + str(Left_Sum) + "\t" \
			+ str(Right_Sum) + "\t" + str(SumDiff) + "\t" + Species_HASH[x] + "\n")

	'''
	Designed for other stats but this should be fine
	'''

	#Get the edge length that proceeds all nodes
	preorder2(n,Association_HASH)
	final_tips = []
	get_tips(n,final_tips)
	print "Tips\t" + str(len(final_tips))
	outw2.write("Tips\t" + str(len(final_tips))+"\n")

	dup = ""
	conf = ""
	conc = ""
	un = ""
	
	for x in Association_HASH:
		if x == "D":
			for i in Association_HASH[x]:
				dup += "\t" + str(i)
			print "Duplication\t" + str(len(Association_HASH[x]))
			outw2.write("Duplication\t" + str(len(Association_HASH[x]))+"\n")
		elif x == "X":
			for i in Association_HASH[x]:
				conf += "\t" + str(i)
			print "Conflict\t" + str(len(Association_HASH[x]))
			outw2.write("Conflict\t" + str(len(Association_HASH[x]))+"\n")
		elif x == "*":
			for i in Association_HASH[x]:
				conc += "\t" + str(i)
			print "Concordance\t" + str(len(Association_HASH[x]))
			outw2.write("Concordance\t" + str(len(Association_HASH[x]))+"\n")
		elif x == "U":
			for i in Association_HASH[x]:
				un += "\t" + str(i)
			print "Uninformative\t" + str(len(Association_HASH[x]))
			outw2.write( "Uninformative\t" + str(len(Association_HASH[x]))+"\n")
