"""This module contains functions for producing a useful, logical output
from the lists of relationships and biparts produced by the other parts
of the program.
"""
import sys

from objects import Node, Rel, Bipart
import comparisons
import read_trees


def sort_conflicts(conflicts):
	"""This should sort a list of biparts that conflict with the species 
	tree into nodes and then from there into separate "conflict topologies".
	It returns a dictionary (complete_dict)  which contains another 
	dictionary for each node, which in turn contains the different 
	conflicts.
	"""

	node_dict = {}
	
	# Add all the conflicts for each node (on the species tree) to a
	# dictionary, sorting them by node id.
	for conflict in conflicts:	
		node = conflict.species_node
		bipart = conflict.ortholog_bipart
		
		if node not in node_dict.keys():
			node_dict[node] = [bipart]
		else:
			node_dict[node].append(bipart)

	# Now we sort each bipart by type of conflict within its node.
	complete_dict = {}

	for key in node_dict.keys():
		bipart_list = node_dict[key]
		conflict_dict = {}
		conflict_dict["conflict_1"] = [bipart_list[0]]

		for bp1 in bipart_list[1:]:
			index = 1	
			
			while True:
				name = "conflict_" + str(index)

				if name in conflict_dict.keys():
					bp2 = conflict_dict[name][0]
				else:
					sys.stderr.write("Making list for node_" + key + " conflict_" + str(index) + "\r")
					conflict_dict[name] = [bp1]
					break

				rel = comparisons.bipart_relationship(bp1, bp2)
				
				if rel == 'concordant':
					conflict_dict[name].append(bp1)
					index = 1
					break
				
				else:
					index += 1

		complete_dict[key] = conflict_dict

	return complete_dict

def length_of_2nd_entry(a_list):
	"""This function takes a list and returns the length of the 2nd entry,
	as long as it is also a list.
	"""
	if isinstance(a_list[1], list):
		return len(a_list[1])
	else:
		print "not a list!"

def conflict_stats(conflicts_dict, tree, outfile):
	"""This function should take a dictionary from sort_conflicts and 
	calculate the most common conflict at each node, second-most common, 
	etc.
	"""
	
	# We made this as a dictionary earlier because it was easier to do it 
	# that way then, but now we want to put things in a defined order so we 
	# need a list. 
	stats_dict = {}

	for node in conflicts_dict.keys():
		stats_dict[node] = []
		
		for conflict in conflicts_dict[node].keys():
			bipart_list = conflicts_dict[node][conflict]
			new_list = [conflict, bipart_list]
			stats_dict[node].append(new_list)
		
	outfile.write("node_id,species_bipart,ortholog_bipart,number_of_conflicts,percentage\n")

	for node in stats_dict.keys():	
		# Order all the conflicts within each node from most to least 
		# common.
		node_on_tree = read_trees.node_finder(tree, node)
		node_bipart = read_trees.postorder3(node_on_tree)
		stats_dict[node].sort(reverse=True, key=length_of_2nd_entry)	

		# Get the total so we can calculate percentages.
		total = 0
		for conflict in stats_dict[node]:
			total += len(conflict[1])

		counter = 0
		cumulative_percent = 0

		for conflict in stats_dict[node]:
			how_common = len(conflict[1])
			percent = float(how_common)/total * 100
			
			# Write each result out to a table. Double-check this!
			output = []
			output.append(str(node))
			output.append(";".join(node_bipart.bipart_proper))
			output.append(";".join(conflict[1][0]))
			output.append(str(how_common))
			output.append(str(percent))
			string = ",".join(output) + "\n"
			outfile.write(string)

			percent = round(percent, 2)
			cumulative_percent += percent
			counter += 1


def make_string(bipart):
	"""This is to make it easier to understand the output when printing 
	things. It makes a neat string out of a bipart made in postorder2 
	that shows both sides.
	"""
	
	name = ''
	bp = bipart.bipart_proper
	not_bp = bipart.other_side

	for i in bp:
		name += i + " "
	
	name += '| '
	
	for i in not_bp:
		name += i + " "

	return str(name)

