#Contains the pxlstr functions 
#Repeats the traversal a lot so possibly could look into fixing that
import numpy as np

#test
def back_home(tree,array):
	
	if tree.parent != None:
		array.append(tree.length)
		back_home(tree.parent,array)




#postorder traverse to get root-to-tip dist
def path_to_tree(tree,lengths):
	
	for child in tree.children:
		
		#get to a child so go back to parent to get val
		if child.istip:

			array = []
			back_home(child,array)
			lengths.append(sum(array))

		path_to_tree(child,lengths)
	
#calculate the variance from a list
#thanks stackoverflow: 35583302
def variance_from_list(results):
	
	#get the mean 
	m = sum(results) / len(results)
	# calculate variance using a list comprehension	
	var_res = sum((xi - m) ** 2 for xi in results) / len(results)
	return var_res

#Use a postorder traversal to get path to tree
def Get_var(Tree,names):
	
	lengths = []
	path_to_tree(Tree,lengths)
	variance = variance_from_list(lengths)
	return variance

def tree_length(tree,array):

	for child in tree.children:
		
		array.append(child.length)
		tree_length(child,array)

	

#call on function to get the tree length
def get_len(gene_root):

	array = []
	tree_length(gene_root,array)
	return sum(array)


'''
#remove the tip on a tree
def remove_tip(tree,taxa):
	
	for child in tree.children:
		
		if child.label == taxa:
			tree.remove_child(child)
		remove_tip(child,taxa)

#If a tip is not a tip then remove
def smoother(tree):

	for child in tree.children:
		
		if len(child.children) == 0:
			
			if child.istip == False:
				tree.remove_child(child.parent.children[0])
				
		smoother(child)

def remove_ogs(tree,array):
	
	for taxa in array:
		remove_tip(tree,taxa)
	smoother(tree)
	if len(tree.children) == 1:
		return tree.children[0]
	else:
		return tree

'''

















	
	
