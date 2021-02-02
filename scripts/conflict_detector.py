import sys
import os
import argparse
import copy

import comparisons
import analysis
import make_trees
import read_trees
from objects import Node, Rel, Bipart

if __name__ == "__main__":
    # Taking user input.
    parser = argparse.ArgumentParser(
        description="A program to help find and analyse conflict and concordance in ortholog trees with duplications.")
    parser.add_argument("-d", "--description", action="store_true",
                        help="Gives a longer description of the program and different modes.")
    parser.add_argument("--species_tree", type=str,
                        help="A file containing the species tree.")
    parser.add_argument(
        "--mode", type=str, help="The mode that we run the program in. Should be n (normal), r (reverse), s (search) or summarize (for analying reverse).")
    parser.add_argument("--gene_tree", type=str,
                        help="File containing a gene tree. Required in 'r' mode.")
    parser.add_argument("--gene_folder", type=str,
                        help="A folder containing all the gene trees to be analysed. Required for 'n' and 's' modes.")
    parser.add_argument("--cutoff", type=int, default=0,
                        help="Nodes with a bootstrap label below the cutoff value will not be included in the analysis.")
    parser.add_argument("--query_bipart", type=str,
                        help="A query to be searched. Required in 's' mode, should be entered in the form: 'taxon1,taxon2', NO SPACES.")
    parser.add_argument("--query_remove", type=str,
                        help="Comma separated taxa you would be ok with not being in the query bipart (only 's' mode). This does not account for tree structure!")
    parser.add_argument("--no_csv", action="store_true",
                        help="Program does not create a .csv file with a more detailed breakdown of the different conflict abundences.")
    parser.add_argument("--outfile_prefix", type=str, default="out", help="Outfile name.")
    parser.add_argument("--output_subtrees", action="store_true", help="Will output subtrees analyzed to a file that has the gene name followed by .subtree")
    parser.add_argument("--annotated_tree", type=str, help="The output of the reverse detector")


    if len(sys.argv) == 1:
        parser.print_usage()
        parser.exit()

    args = parser.parse_args()
    if args.description:
        print "This program shows gene tree conflict whilst accounting\n\
for duplication events in gene trees.\n\n\
Normal mode (specify mode as 'n') produces two outputs:\n\
the printed output shows two trees in newick format,\n\
showing the number of conflicts and concordances at each\n\
node as labels. It also outputs a tab-separated table\n\
showing a more detailed breakdown of the types of\n\
conflict at each node.\n\n\
reverse mode (specify mode as 'r') maps conflicts with\n\
the species tree back onto the user-specified gene tree,\n\
labelling concordance with '*', conflict with 'X' and a\n\
duplication node with 'D'.\n\n\
Search mode (specify mode as 's') allows the user to\n\
search for a specific bipart in a folder of gene trees,\n\
and returns a list of all the tree files that contain\n\
that bipart.\n\nmode summarize takes in the annotated tree\n\
from -r mode and will summarize the results."
        sys.exit(0)

    mode = args.mode
    species_tree = args.species_tree
    cutoff = args.cutoff
    gene_folder = args.gene_folder
    gene_tree = args.gene_tree
    query_bipart = args.query_bipart
    outfile_prefix = args.outfile_prefix
    outfile_subtrees = args.output_subtrees
    query_remove = args.query_remove
    annotated_tree = args.annotated_tree

    if mode == 'n':
        # Check we have arguments required.
        if not species_tree or not gene_folder:
            print "--species_tree and --gene_folder are required arguments in this mode."
            sys.exit(0)

        # Making the species tree.
        tree_file = open(species_tree, "r")
        for line in tree_file:
            tree = line
        species_root, species_name_array = make_trees.build(tree)
        species_biparts = read_trees.postorder2(species_root)
        all_taxa = read_trees.postorder3(species_root)
        species_biparts.append(all_taxa)

        # We need these later.
        total_conflicts = []
        total_concordances = []

        # Making sure folder always has the same name.
        if gene_folder[-1] == '/':
            homologs_folder = gene_folder
        else:
            homologs_folder = gene_folder + '/'

        # Because we do all the work inside this for loop, we use less
        # memory as we only handle one file at once.
        file_list = os.listdir(homologs_folder)
        len_file_list = str(len(file_list))
        file_no = 0
        for file in file_list:
            file_no += 1
            sys.stderr.write("Processing file " + str(file_no) +
                             " of " + len_file_list + ".\r")

            # Building the gene tree.
            file_location = str(homologs_folder) + str(file)
            gene_file = open(file_location, "r")


            subtree_printout = []
            printout_filename = []
            for line in gene_file:
                tree = line
                gene_root, gene_name_array = make_trees.build(tree)

                # We make 'subtrees' to deal with the problem of
                # duplication: the subtrees contain no duplications and
                # cover the whole gene tree.
                trees = make_trees.subtrees_function(gene_root)
                log_name = gene_folder[:-1] + "_subtree_comp.log"

                # This will act on just orthologs.
                if len(trees) == 1:
                    for tree in trees:
                        conflicts, concordances = comparisons.compare_trees(
                            species_biparts, species_name_array, tree, [], "ortholog", mode, "some_log_name", cutoff, str(file))
                        total_concordances.extend(concordances)

                        # We only use one conflict per gene (accounting for nesting).
                        filtered_conflicts = comparisons.filter_conflicts(
                            conflicts)
                        total_conflicts.extend(filtered_conflicts)
                else:
                	'''
                	Procedure 1) Make Tree 2) Identify duplication 3) Cut out subtrees
                	based on duplications 4) For main tree and subtrees collapse into
                	polytomies with a note not to analyze them 5) Analyze the tree set
                	'''
                	#1) Make Tree
                	g, g_name = make_trees.build(tree)
                	
                	#2) Identify the duplications
                	make_trees.label_duplications(g)

                	# make a tree array that will contain all subtrees divided out
                	# based on dups and polytomies unanalyzed because of downstream
                	# procedure
                	trarray = []
                	#extra names should be the ones applicable to all subtrees
                	extra_names = []
                	#3) Cut out subtrees based on duplications (e.g left and right of dup)
                	make_trees.subtree_divide(g, trarray, extra_names)                	
                	trarray.append(g)
                	extra_names.append([]) #final tree has all the names anyway

                	#Loop to process all subtrees and the main tree
                	count = -1
                	for tree in trarray:
						
						count += 1
						#4) collapse duplications into polytomies
						dup =[]
						dup = tree.get_dup_count()
						#this assumes no dups overlap so it's probably overkill but
						#should work
						for i in dup:	
							tree.collapse_dups()
						
						if tree.istip == True:
							pass
						#this gets passed because these are counted as part of the larger tree
						elif tree.label == "D":
							pass
						#5) Analyze the subtree
						else:

							conflicts, concordances=comparisons.compare_trees(
							    species_biparts, species_name_array, tree, extra_names[count], "homolog", mode, "some_log_name", cutoff, str(file))
							total_concordances.extend(concordances)

							filtered_conflicts=comparisons.filter_conflicts(conflicts)
							total_conflicts.extend(filtered_conflicts)
							if outfile_subtrees:
								tree.get_rid_of_note("CollapsedNotCounted")
								#this is actually taxa on the other side of clade
								printout_filename.append(extra_names[count])
								subtree_printout.append(tree)
			if outfile_subtrees:
				outw = open(file+".subtree","w")
				for x in range(0,len(subtree_printout)):
					subtree_printout[x].get_rid_of_note("CollapsedNotCounted")
					make_trees.add_loci(subtree_printout[x])
					if len(printout_filename[x]) != 0:
						t = "((" + ",".join(printout_filename[x]) + ")-1," + subtree_printout[x].get_newick_repr() + ");"
					else:
						t = subtree_printout[x].get_newick_repr() + ";"
					outw.write(t+"\n")
				outw.close()					

        if not args.no_csv:
            # Extra analysis to get the relative frequenices of the
            # conflicts, etc.
            outfile=open(str(outfile_prefix) + "_analysis.csv", "w")
            sorted_conflicts=analysis.sort_conflicts(total_conflicts)
            analysis.conflict_stats(sorted_conflicts, species_root, outfile)

        concord_out = open(str(outfile_prefix) + "_concord.tre", "w")
        conflict_out = open(str(outfile_prefix) + "_conflict.tre", "w")
        labels_out = open(str(outfile_prefix) + "_labels.tre", "w")
        concord_gene_out = open(str(outfile_prefix) + "_concord_gene_counts.tre", "w")
        conflict_gene_out = open(str(outfile_prefix) + "_conflict_gene_counts.tre", "w")
        if cutoff == 0:
        	nodes_analyzable_out = open(str(outfile_prefix) + "_total_analyzed.tre", "w")
        

        # Map concordances and conflicts back onto the tree using the lists
        print "\n"
        make_trees.clear_labels(species_root)
        make_trees.tree_map(species_root, total_concordances)        
        make_trees.add_zeros(species_root)
        concordance_tree=species_root.get_newick_repr(showbl=True)
        print "Concordance tree: "
        print concordance_tree + ";"
        concord_out.write(concordance_tree + ";")		

        make_trees.clear_labels(species_root)
        make_trees.tree_map(species_root, total_conflicts)
        make_trees.add_zeros(species_root)
        conflict_tree=species_root.get_newick_repr(showbl=True)
        print "Conflict tree:"
        print conflict_tree + ";"
        conflict_out.write(conflict_tree + ";")
        
        if cutoff == 0:
            total_rels = []
            total_rels.extend(total_concordances)
            total_rels.extend(total_conflicts)
            make_trees.clear_labels(species_root)
            make_trees.tree_map(species_root, total_rels)
            make_trees.add_zeros(species_root)
            analyzed_tree = species_root.get_newick_repr(showbl=True)
            nodes_analyzable_out.write(analyzed_tree + ";")
        	

        make_trees.label_node_ids(species_root)
        labels_tree=species_root.get_newick_repr(showbl=True)
        print "Labels tree: "
        print labels_tree + ";"
        labels_out.write(labels_tree + ";")
    
    	#As opposed to all the concordances this only analyzes genes
        concord_gene_list = str(outfile_prefix) + "_concord_genes.csv"
        concord_gene_dup_out = str(outfile_prefix) + "_multi_concord_gene_counts.csv"
        make_trees.clear_labels(species_root)
        make_trees.tree_map3(species_root, total_concordances, concord_gene_list, concord_gene_dup_out)
        make_trees.add_zeros(species_root)
        concordance_gene_tree=species_root.get_newick_repr(showbl=True)
        print "Concordant gene tree: "
        print concordance_gene_tree + ";"
        concord_gene_out.write(concordance_gene_tree + ";" + "\n")
        concord_gene_out.close()
        
        #As opposed to all the conflicts this only analyzes genes
        conflict_gene_list = str(outfile_prefix) + "_conflict_genes.csv"
        conflict_gene_dup_out = str(outfile_prefix) + "_multi_conflict_gene_counts.csv"
        make_trees.clear_labels(species_root)
        make_trees.tree_map3(species_root, total_conflicts, conflict_gene_list, conflict_gene_dup_out)
        make_trees.add_zeros(species_root)
        conflict_gene_tree=species_root.get_newick_repr(showbl=True)
        print "Conflict gene tree: "
        print conflict_gene_tree + ";"
        conflict_gene_out.write(conflict_gene_tree + ";" + "\n")
        conflict_gene_out.close()
        
    #This is to map species tree to gene tree
    elif mode == "r":
		
		#check for args
        if not species_tree or not gene_tree:
            print "--species_tree and --gene_tree are required arguments in this mode."
            sys.exit()
        
        #read in the species tree
        file = open(species_tree, "r")
        for line in file:
        	tree = line
        
        #get species tree as node, get biparts of species tree
        #get bipart for root
        species_root, species_name_array = make_trees.build(tree)
        species_biparts = read_trees.postorder2(species_root)
        all_taxa = read_trees.postorder3(species_root)
        species_biparts.append(all_taxa)
        
        #read in gene tree and make sure only one exists
        file = open(gene_tree, "r")
        tree = file.readline()
    	if file.readline():
    		print "There must be only one tree in the gene tree file!"
    		sys.exit()
    	
    	#turn gene tree into node structure
    	tree_root, gene_name_array = make_trees.build(tree)
    	g, g_array = make_trees.build(tree)
    	
    	'''
    	1) Map duplication nodes 2) Split into subtrees based on duplication nodes
    	3) Analyze the subtrees for conflict and concordance 4) Blank nodes are changed
    	to uninformative
    	'''
    	
    	#1) Map duplication nodes
    	make_trees.label_duplications(tree_root)
    	
    	#first tree gets mutated and loses node positions during 
    	#the subtree procedure so we need a second to map back onto
    	make_trees.label_duplications(g)
    	
    	#2) divide into subtrees no deep copy so any manipulation will affect it
    	trarray = []
    	extra_names = []
    	make_trees.subtree_divide(tree_root,trarray,extra_names)
    	trarray.append(tree_root)
    	extra_names.append([])
    	#print extra_names
    	if outfile_subtrees:
    		tree_to_print = []
    		out_subs = open(gene_tree+".subtree","w")
    	
    	#parse through the subtrees
    	count = -1
    	for tree in trarray:
			count += 1
			#A slight complexity is that it will be collapsed into subtrees
			#This will lose node ID and total taxa, so this takes in total taxa
			pre_col_taxa = tree.lvsnms()
			
    		#collapse duplications into polytomies, these are the trees
    		#with necessary sampling to check for conflict and concordance
			dup =[]
			dup = tree.get_dup_count()
			for i in dup:	
				tree.collapse_dups()
			
			#tips cannot be labeled
			if tree.istip == True:
				pass
			
			#duplications cannot be labeled
			elif tree.label == "D":
				pass

			#here are the trees that need to be analyzed against the species tree
			else:

				#get all taxa in subtree for the rel analysis
				subtree_taxa = []
				subtree_taxa = tree.lvsnms()
				#print "Tree is: " + tree.get_newick_repr()
				#print "Taxa in tree: " + str(subtree_taxa)
				#print extra_names[count]
				subtree_taxa.extend(extra_names[count])
				
				#print "Added taxa: " + str(subtree_taxa)
				
				#get the biparts
				subtree_biparts = read_trees.postorder2(tree)
				all_taxa = read_trees.postorder3(tree)
				subtree_biparts.append(all_taxa)
				
				#get how things are related to each other
				rel_list = []
				
				#Uniformative is chosen a different way here so always give it no cutoff
				rel_list = comparisons.comp_biparts(subtree_biparts, species_biparts, subtree_taxa, species_name_array, "some_log_name", 0, "r", "")

				#divide up conflicts and concordances on the tree
				conflicts = []
				concordances = []

				#This identifies conflicts and concordances on the subtrees
				for rel in rel_list:
					if rel.relation == 'conflict':
						conflicts.append(rel)
					elif rel.relation == 'concordant':
						concordances.append(rel)
				
				#find the subtree we are at in the main tree, mrca will become a mutable
				#part of the main tree so it should co-label as subtrees are analyzed
				mrca = Node()
				
				#find the clade this exists in
				g.find_mrca(pre_col_taxa,mrca)
			
				#4) Identify what is a concordance, conflict, uninformative 
				#and label using the label node function within the Node structure
				for x in concordances:
					mrca.label_node(x.ortholog_bipart,"*",cutoff)
				for x in conflicts:
					mrca.label_node(x.ortholog_bipart,"X",cutoff)
				
				#for printing subtrees
				if outfile_subtrees:
					out_mrca = copy.deepcopy(tree)
					out_mrca.clr_label()
					make_trees.add_loci(out_mrca)
					if len(extra_names[count]) != 0:
						out_subs.write("((" + ",".join(extra_names[count]) + ")," + out_mrca.get_newick_repr() + ");\n")
					else:
						out_subs.write(out_mrca.get_newick_repr() + ";\n")
        #add the loci ID's back on
        make_trees.add_loci(g)
        outw = open(str(outfile_prefix) + "_concon.tre", "w")
        outw.write(g.get_newick_repr(showbl=True)+";")
        outw.close()
					
    #Mode to search for a bipart
    #Option 1 allow missing taxa assuming they are missing in the subtree
    #Option 2 find only subtrees that contain all taxa for a bipart
    #Currently neither is implemented
    elif mode == "s":
		
		#test bipart MJM2704_Achatocarpus_gracilis,MJM1677_Phaulothamnus_spinescens
		# Checking arguments.
		if not gene_folder or not query_bipart:
			print "--gene_folder and --query_bipart are required arguments in this mode."
			sys.exit()

		# We need these for later.
		concordant_files=[]

		# Putting the query into a form comp_biparts likes.
		query_entered = query_bipart.split(",")

		# Opening and reading the provided folder of homologs.
		if gene_folder[-1] == "/":
			homologs_folder = gene_folder
		else:
			homologs_folder = gene_folder + "/"
		
		file_list = os.listdir(homologs_folder)

		count = 0
		subtree_printout = []
		printout_filename = []
		final_print_matches = []
		final_print_names = []
		
		if query_remove:
			
			removable = query_remove.split(",")

		#This becomes a basic subtree search method
		if outfile_subtrees:
				outw = open(outfile_prefix+".search.subtree","w")
		printout_filename = []
		subtree_printout = []
		names_on_otherside = []
		
		
		for filename in file_list:
			
			count += 1
			
			total_concordances = []
			sys.stderr.write("processing " + str(count) + " of " + str(len(file_list)) + '\r')

			# Building a tree and making subtrees.
			file_location = str(homologs_folder) + str(filename)
			file = open(file_location, "r")
            
			tree = file.readline()
			if file.readline():
				print "WARNING!!!! This is only analyzing the first tree in your file"
			
			#build the tree
			gene_root, gene_name_array = make_trees.build(tree)
			
			'''
			This is only concordance which makes life a lot easier
			so 1) map the duplications 2) divide into subtrees based on 
			duplications 3) send subtrees and their total taxa into a 
			concordance identified
			'''
			#1) label the duplication
			make_trees.label_duplications(gene_root)
			
			#2) divide into subtrees
			trarray = []
			extra_names = []
			make_trees.subtree_divide(gene_root,trarray,extra_names)
			trarray.append(gene_root)
			#the final item in extra names is going to be everything that has been collapsed
			root_names = []
			list_set = set()
			for x in extra_names:
				for y in x:
					list_set.add(y)
			all_names = gene_root.lvsnms()
			root_names = (set(all_names) - list_set)
					
			
			extra_names.append(list(root_names))
			
			#this will have the biparts that are found
			ret_bipart = []
			count2 = -1
			
			for tree in trarray:
				
				count2 += 1
				dup =[]
				dup = tree.get_dup_count()
				
				for i in dup:	
					tree.collapse_dups()

        		#tips cannot be labeled
				if tree.istip == True:
					pass
			
				#duplications cannot be labeled
				elif tree.label == "D":
					pass

				#here are the trees that need to be analyzed against the species tree
				else:
				
					'''
					Missing taxa search algorithm. 3a) Take the taxa that can be removed and check if 
					they are all in the subtree. 3b) If they are then great and move on since an exact 
					match could exist. If they are not then take the ones missing out of the query.
					3c) If the query still has at least two taxa (allowing the rel to exist) then proceed
					with just the two as the query bipart. if less than 2 exist then don't proceed
					3d) Since extra names is essentially the outgroup bipart then it should check to make
					sure the bipart is not split up if you've allowed missing taxa
					'''
					search_list = []
					#do non-exact search
					if query_remove:
						
						#3a) check if removable taxa are in subtree
						
						#get all taxa in subtree
						subtree_taxa = tree.lvsnms()
						
						#compare to taxa that can be removed, not_in_subtree will have the
						#ones that can be removed
						#3d) make sure no intersection exists between the outgroup and the query
						check =  any(item in query_entered for item in extra_names[count2])
						if check:
							pass
						else:
							not_in_removable = []
							not_in_subtree = []
							not_in_removable,not_in_subtree = comparisons.unique_array(removable,subtree_taxa)
						
							#print "Subtree list: " + str(subtree_taxa) + "\tRemovable\t" + str(removable) 
							#print not_in_subtree
						
							#3b) remove the ones not in the subtree from the original query to get a new
							#search list
							search_list = [x for x in query_entered if x not in not_in_subtree]
						
							#print "For searching: " + str(search_list)
						
						
					else:
					
						search_list = query_entered

					
					#find concordance, give a subtree, the total taxa in the subtree
					#the query

					cur_len = len(ret_bipart)
					#Make sure you're searching for a bipart and not a tip or nothing
					if 1 < len(search_list):
						tree.find_bipart(search_list,cutoff,ret_bipart)
					
					#If the length of biparts has expanded print out the subtree
					if cur_len != len(ret_bipart):
						if outfile_subtrees:
							tree.get_rid_of_note("CollapsedNotCounted")
							printout_filename.append(filename)
							out_mrca = copy.deepcopy(tree)
							make_trees.add_loci(out_mrca)
							subtree_printout.append(out_mrca)
							names_on_otherside.append(extra_names[count2])
						
					
			if len(ret_bipart) != 0:
			#	print filename + " has " + str(len(ret_bipart))
				final_print_names.append(filename)
				final_print_matches.append(str(len(ret_bipart)))
					
        
			#if count == 10:
			#	sys.exit()
		if outfile_subtrees:
			for x in range(0,len(subtree_printout)):
				tree.get_rid_of_note("CollapsedNotCounted")
				make_trees.add_loci(tree)
				#t = printout_filename[x] + ": ((" + subtree_printout[x].get_newick_repr(showbl = True) + "," + ",".join(names_on_otherside[x]) + "));"
				t = printout_filename[x] + ": " + subtree_printout[x].get_newick_repr(showbl = True)
				outw.write(t+"\n")
			outw.close()
    	
		total_hits = 0
		outw = open(outfile_prefix + "_search.csv", "w")
		for x in range(0,len(final_print_names)):
			outw.write(str(final_print_names[x]) + "," + str(final_print_matches[x])+"\n")
			total_hits += int(final_print_matches[x])
		outw.close()
		print "total gene families: " + str(len(final_print_names)) + " total matches: " + str(total_hits) 
   
    elif mode == "summarize":
		if not annotated_tree:
			print "--annotated_tree is required in this mode."
			sys.exit(0)
		file = open(annotated_tree, "r")
		tree = file.readline()
		gene_root, gene_name_array = make_trees.build(tree)
		counts = []
		possible = ["Duplications","D","Concordances","*","Conflicts","X","Uninformative","U"]
		gene_root.count_label(counts)
		for i in range(0,len(possible),2):
			print possible[i] + ": " + str(counts.count(possible[i+1]))
		print "Total tips: " + str(len(gene_name_array))
			
		if query_bipart and species_tree:
			if outfile_prefix:
				outfile_sum = open(outfile_prefix + ".tsv", "w")
			else:
				outfile_sum = open("summary.tsv", "w")
			
			'''
			Gonna be a pain but same ideas as above. Tree should already have duplications labeled
			so its a matter of identifying the subtrees and no collapsing said duplications
			1) Divide species tree into biparts 2) Divide at the duplication into subtrees 
			3) get the stats on child 1 and child 2 of each subtree 4) profit
			The key is the ingroup gives the point of duplication on the species tree (e.g all unique taxa)
			The outgroup determines if said duplication is concordant with the species tree
			'''
			file2 = open(species_tree,"r")
			sp_tr = file2.readline()
			#1) divide species tree into biparts
			bipart_array = []
			species_root, species_names = make_trees.build(sp_tr)
			species_root.get_biparts(bipart_array)
			
			query_entered = query_bipart.split(",")
			
			#2) divide duplications into subtrees
			trarray = []
			extra_names = []
			make_trees.subtree_divide_at_base_of_dup(gene_root,trarray,extra_names)
			
			count2 = -1
			outfile_sum.write("LeftTips\tLeftDuplications\tLeftConcord\tLeftConflict\tLeftUninformative\tRightTips\tRightDuplications\tRightConcord\tRightConflict\tRightUninformative\n")
			for tree in trarray:
				count2 += 1
				
				if tree.istip == False and tree.label == "D":
					#print "Here is the tree: " + str(tree.lvsnms_uniq())
					#print "Here is the names: " + str(extra_names[count2])
					
					#First need to check that the query bipart isn't in the extra names
					#e.g. the query doesn't intersect with outgroups creating conflict
					check =  any(item in query_entered for item in extra_names[count2])
					if check:
						pass
					else:
						ingroup = tree.lvsnms_uniq()
						#print "Ingroup: " + str(ingroup) + " Outgroup: " + str(extra_names[count2])
						#Get the species tree node the ingroup is associated with (e.g smallest
						#bipart to have all taxa in it)
						sp_dup_node = []
						sp_dup_node = comparisons.get_smallest_match(bipart_array,ingroup)
						#print "Ingroup: " + str(ingroup) + " duplicated at " + str(sp_dup_node)
						
						#The species node matches the query entered, here the analysis can
						#be conducted
						if sorted(query_entered) == sorted(sp_dup_node):
							child1_counts = []
							child2_counts = []
							tree.children[0].count_label(child1_counts)
							tree.children[1].count_label(child2_counts)
							#print "Left side: " + str(len(tree.children[0].lvsnms())) + " Right side: " + str(len(tree.children[1].lvsnms()))
							#print str(child1_counts) + " compared to " + str(child2_counts)
							outfile_sum.write(str(len(tree.children[0].lvsnms())) + "\t" + str(child1_counts.count(possible[1])) \
								+ "\t" + str(child1_counts.count(possible[3])) + "\t" + str(child1_counts.count(possible[5])) \
								+ "\t" + str(child1_counts.count(possible[7])) + "\t" + str(len(tree.children[1].lvsnms())) \
								+ "\t" + str(child2_counts.count(possible[1])) + "\t" + str(child2_counts.count(possible[3])) \
								+ "\t" + str(child2_counts.count(possible[5])) + "\t" + str(child2_counts.count(possible[7])) + "\n") 
					
		else:
			print "For more detailed analysis give species tree and a bipartition of interest"			
			
		
    	
    
    else:
    	print "Mode not found, for options type: python Conflict_detector.py -h"   
   
    		
    '''
    elif mode == "r":
        # Checking we have the right arguments.
        if not species_tree or not gene_tree:
            print "--species_tree and --gene_tree are required arguments in this mode."
            sys.exit()

        # Again, we start by reading in the species tree.
        file=open(species_tree, "r")
        for line in file:
            tree=line
        species_root, species_name_array=make_trees.build(tree)
        species_biparts=read_trees.postorder2(species_root)
        all_taxa=read_trees.postorder3(species_root)
        species_biparts.append(all_taxa)

        # Then we read in the gene tree.
        file=open(gene_tree, "r")
        counter=0
        for line in file:
            tree=line
            counter += 1
        if counter != 1:
            print "There must be only one tree in the gene tree file!"
            sys.exit()
        gene_root, gene_name_array=make_trees.build(tree)
        trees=make_trees.subtrees_function(gene_root)

        if len(trees) != 0:
            outfile=gene_tree + "_subtree_comp.log"
            for tree in trees:
                conflicts, concordances=comparisons.compare_trees(
                    species_biparts, species_name_array, tree, mode, outfile, cutoff)
                make_trees.tree_map2(tree, conflicts, 'X')
                make_trees.tree_map2(tree, concordances, '*')

            # As in normal mode, there are some nodes missed by the
            # subtree method.
            tricky_nodes=read_trees.identify_tricky_nodes(gene_root, trees)
            conflicts=[]
            concordances=[]

            for node in tricky_nodes:
                node_bipart=read_trees.postorder3(node)
                node_bipart_list=[node_bipart]
                node_name_array=read_trees.postorder3(node)
                node_name_array=node_name_array.bipart_proper
                new_names=[]
                current_node=node

                while True:
                    parent=current_node.parent
                    if parent != None:
                        make_trees.label_duplications(parent, recursive=False)
                        if parent.label == 'D':
                            current_node=parent
                        elif type(parent) is None:
                            break
                        else:
                            for child in parent.children:
                                if child != current_node:
                                    new_names=read_trees.postorder3(child)
                                    new_names=new_names.bipart_proper
                            break
                    else:
                        break

                node_name_array.extend(new_names)
                outfile=gene_tree + "_tricky_nodes_comp.log"
                rel_list=comparisons.comp_biparts(
                    node_bipart_list, species_biparts, node_name_array, species_name_array, outfile, cutoff, mode)

                for rel in rel_list:
                    if rel.relation == 'conflict':
                        conflicts.append(rel)
                    elif rel.relation == 'concordant':
                        concordances.append(rel)

            # Printing the output.
            make_trees.tree_map2(gene_root, conflicts, 'X')
            make_trees.tree_map2(gene_root, concordances, '*')
            make_trees.label_duplications(gene_root)
            make_trees.label_uninformative(gene_root)
            make_trees.add_loci(gene_root)
            new_tree=gene_root.get_newick_repr(showbl=True)
            # print new_tree + ";"
            if outfile_prefix:
                out=open(str(outfile_prefix) + "_concon.tre", "w")
                labels=open(str(outfile_prefix) + "_labels.tre", "w")
            else:
                out=open("concon.tre", "w")
                labels=open("labels.tre", "w")
            out.write(new_tree + ";")
            make_trees.label_node_ids(gene_root)
            labels_tree=gene_root.get_newick_repr(showbl=True)
            labels.write(labels_tree + ";")
        else:
            # It's a lot simpler when there are no duplications.
            outfile=gene_tree + "_comp.log"
            conflicts, concordances=comparisons.compare_trees(
                species_biparts, species_name_array, gene_root, mode, outfile, cutoff, mode)
            make_trees.tree_map2(gene_root, conflicts, 'X')
            make_trees.tree_map2(gene_root, concordances, '*')
            make_trees.label_duplications(gene_root)
            make_trees.label_uninformative(gene_root)
            make_trees.add_loci(gene_root)
            new_tree=gene_root.get_newick_repr(showbl=True)
            print new_tree + ";"
            if outfile_prefix:
                out=open(str(outfile_prefix) + "_concon.tre", "w")
                labels=open(str(outfile_prefix) + "_labels.tre", "w")
            else:
                out=open("concon.tre", "w")
                labels=open("labels.tre", "w")
            out.write(new_tree + ";")
            make_trees.label_node_ids(gene_root)
            labels_tree=gene_root.get_newick_repr(showbl=True)
            labels.write(labels_tree + ";")

    elif mode == 's':
        # Checking arguments.
        if not gene_folder or not query_bipart:
            print "--gene_folder and --query_bipart are required arguments in this mode."
            sys.exit()

        # We need these for later.
        concordant_files=[]

        # Putting the query into a form comp_biparts likes.
        query_entered=query_bipart.split(",")
        query=Bipart("query", "query")
        for taxon in query_entered:
            query.add_taxon(taxon)
        query_name_array=query_entered

        # Opening and reading the provided folder of homologs.
        if gene_folder[-1] == "/":
            homologs_folder=gene_folder
        else:
            homologs_folder=gene_folder + "/"
        file_list=os.listdir(homologs_folder)

        for filename in file_list:
            total_concordances=[]
            sys.stderr.write("processing " + str(filename) + '\r')

            # Building a tree and making subtrees.
            file_location=str(homologs_folder) + str(filename)
            file=open(file_location, "r")
            for line in file:
                tree=line
            gene_root, gene_name_array=make_trees.build(tree)
            trees=make_trees.subtrees_function(gene_root)

            # Looking for concordances in the subtrees.
            outfile=gene_folder[:-1] + "_comp.log"
            for tree in trees:
                concordances=[]
                conflicts, concordances=comparisons.compare_trees(
                    [query], query_name_array, tree, mode, outfile, cutoff)
                total_concordances.extend(concordances)

            # Dealing with the nodes that aren't in subtrees again.
            tricky_nodes=read_trees.identify_tricky_nodes(gene_root, trees)
            for node in tricky_nodes:
                concordances=[]
                node_bipart=read_trees.postorder3(node)
                node_bipart_list=[node_bipart]
                node_name_array=read_trees.postorder3(node)
                node_name_array=node_name_array.bipart_proper

                new_names=[]
                current_node=node
                while True:
                    parent=current_node.parent
                    if parent == None:
                        break

                    make_trees.label_duplications(parent, recursive=False)
                    if parent.label == 'D':
                        current_node=parent
                    else:
                        for child in parent.children:
                            if child != current_node:
                                new_names=read_trees.postorder3(child)
                                new_names=new_names.bipart_proper
                        break

                node_name_array.extend(new_names)
                outfile=gene_folder[:-1] + "_tricky_nodes_comp.log"
                rel_list=comparisons.comp_biparts(
                    [query], node_bipart_list, query_name_array, node_name_array, outfile, cutoff, mode)

                for rel in rel_list:
                    if rel.relation == 'concordant':
                        concordances.append(rel)

                total_concordances.extend(concordances)

            # We add all files with any concordances to the list.
            for i in range(len(total_concordances)):
                concordant_files.append(filename)

        # Printing the list of files concordant.
        print "\n"
        for file in concordant_files:
            print str(file)

        number_of_times=len(concordant_files)
        concordant_set=set(concordant_files)
        total_files=len(concordant_set)

        # Printing the total length of the list.
        print "Relationship " + str(query_bipart) + " appears " + str(number_of_times) + \
            " times in " + str(total_files) + " genes of " + \
            str(len(file_list)) + " total genes."
    '''