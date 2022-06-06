#!/usr/bin/python
# coding=utf-8
import sys
import os
import argparse
import copy

import comparisons
import analysis
import make_trees
import read_trees
import sortadate
from objects import Node, Rel, Bipart

if __name__ == "__main__":
    # Taking user input.
    parser = argparse.ArgumentParser(
        description="A program to help find and analyse conflict and concordance in ortholog trees with duplications.")
    parser.add_argument("-d", "--description", action="store_true",
                        help="Gives a longer description of the program and different modes.")
    parser.add_argument("--species_tree", type=str,
                        help="A file containing the species tree. Required in 'n' and 'r' modes. All trees should be in newick format.")
    parser.add_argument(
        "--mode", type=str, help="The mode that we run the program in. Should be n (normal), r (reverse), s (search), summarize (for analying reverse) or sortadate (gene filtering, does not remove outgroups).")
    parser.add_argument("--gene_tree", type=str,
                        help="File containing a gene tree. Required in 'r' mode. All trees should be in newick format.")
    parser.add_argument("--gene_folder", type=str,
                        help="A folder containing all the gene trees to be analysed. Required for 'n' and 's' modes. All trees should be in newick format.")
    parser.add_argument("--cutoff", type=int, default=0,
                        help="Nodes with a bootstrap label below the cutoff value will not be included in the analysis.")
    parser.add_argument("--query_bipart", type=str,
                        help="A query to be searched. Required in 's' mode, should be entered in the form: 'taxon1,taxon2', NO SPACES.")
    parser.add_argument("--query_remove", type=str,
                        help="Comma separated taxa you would be ok with not being in the query bipart (only 's' mode). This does not account for tree structure!")
    parser.add_argument("--no_csv", action="store_true",
                        help="Program does not create a .csv file with a more detailed breakdown of the different conflict abundances in normal mode.")
    parser.add_argument("--outfile_prefix", type=str, default="out", help="A prefix for all of the outfiles produced by the program in whichever mode you run it in. ")
    parser.add_argument("--output_subtrees", action="store_true", help="Will output subtrees analyzed to a file that has the gene name followed by .subtree.")
    parser.add_argument("--annotated_tree", type=str, help="File containing a gene tree labelled by 'r' mode. Required in 'summarize' mode.")
    parser.add_argument("--perfect_concordance", action="store_true", help="If this argument is used in normal mode, the program will give a list of all the trees in the input folder that contain 0 conflicts with the species tree.")
    parser.add_argument("--perfect_concordance_strict", action="store_true", help="Same as perfect_concordance, except that only gene trees containing every taxon within the species tree will be included.") 

    if len(sys.argv) == 1:
        parser.print_usage()
        parser.exit()

    args = parser.parse_args()
    if args.description:
        print "This program identifies gene tree conflict whilst accounting\n\
for duplication events in gene trees.\n\n\
Normal mode (specify mode as 'n') maps conflicts and \n\
concordances from a folder of gene trees (one per file)\n\
back onto a species tree, and reports the frequencies of \n\
conflicts and concordances at each node in the species \n\
tree. \n\n\
Reverse mode (specify mode as 'r') maps conflicts with\n\
the species tree back onto the user-specified gene tree,\n\
labelling concordance with '*', conflict with 'X' and a\n\
duplication node with 'D'.\n\n\
Search mode (specify mode as 's') allows the user to\n\
search for a specific bipart in a folder of gene trees,\n\
and returns a list of all the tree files that contain\n\
that bipart.\n\n\
Summarize mode (specify mode as 'summarize') reports basic \n\
statistics for an annotated tree, such as the total number \n\
of duplications, conflicts and concordances for that tree. \n\n\
For more information, including a more detailed explanation \n\
of the output files, please see the manual at [PLACEHOLDER \n\
TEXT]. \n"
        sys.exit(0)

    mode = args.mode
    species_tree = args.species_tree
    cutoff = args.cutoff
    gene_folder = args.gene_folder
    gene_tree = args.gene_tree
    query_bipart = args.query_bipart
    outfile_prefix = args.outfile_prefix
    outfile_subtrees = args.output_subtrees # ???
    query_remove = args.query_remove
    annotated_tree = args.annotated_tree
    
    # If return perfect concordance strict, we still want to return perfect 
    # concord.
    ret_perf_concord = args.perfect_concordance
    ret_perf_concord_strict = args.perfect_concordance_strict
    if ret_perf_concord_strict:
        ret_perf_concord = True

    # Mapping conflicts and concordances onto the species tree. 
    if mode == 'n':
        # Check we have arguments required.
        if not species_tree or not gene_folder:
            print "--species_tree and --gene_folder are required arguments in \
this mode."
            sys.exit(0)

        # In 'return perfect concordances' mode, we need to initialise a list
        # to hold the files
        if ret_perf_concord:
            perf_concords = []

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

        # Making sure the folder name is correct.
        if gene_folder[-1] == '/':
            homologs_folder = gene_folder
        else:
            homologs_folder = gene_folder + '/'

        # Because we do all the work of counting conflicts inside this for loop,
        # we use less memory as we only handle one file at once.
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

            # For later, if outfile_subtrees is input. Gives detailed 
            # information on the subtrees. 
            subtree_printout = []
            printout_filename = []

            # Reading the gene_file. Can be one or many lines?
            for line in gene_file:
                tree = line
                gene_root, gene_name_array = make_trees.build(tree)

                # We make 'subtrees' to deal with the problem of duplication:
                # the subtrees will ultimately contain no duplications and cover
                # the whole gene tree.
                trees = make_trees.subtrees_function(gene_root)

                # This will act on trees with no duplications (orthologs).
                if len(trees) == 1:
                    for tree in trees:
                        conflicts, concordances = comparisons.compare_trees(
                            species_biparts, 
                            species_name_array, 
                            tree, 
                            [], 
                            "ortholog", 
                            mode, 
                            "some_log_name", 
                            cutoff,
                            str(file))
                        total_concordances.extend(concordances)

                        # We only use one conflict per gene (accounting for 
                        # nesting).
                        filtered_conflicts = comparisons.filter_conflicts(
                            conflicts)
                        
                        total_conflicts.extend(filtered_conflicts)
                        
                        # If we're just looking for the perfectly concordant ones
                        if ret_perf_concord:
                            if len(filtered_conflicts) != 0:
                                perfect = False
                            else:
                                perfect = True

                else:
                    # Make a second gene tree (copy).
                    g, g_name = make_trees.build(tree)
                    
                    # Label duplications as the subtrees are based on 
                    # duplications.
                    make_trees.label_duplications(g)
                    
                    trarray = []
                    extra_names = []
                    
                    # The subtrees are "cut out" based on duplications and go 
                    # into trarray/extra_names. 
                    make_trees.subtree_divide(g, trarray, extra_names)                  
                    trarray.append(g)  
                    extra_names.append([]) #final tree has all the names anyway

                    if ret_perf_concord:
                        perfect = True

                    # Processing the subtrees.
                    count = -1
                    for tree in trarray:
                        count += 1

                        # Collapsing the duplications into polytomies means they
                        # aren't analysed in that subtree. 
                        dup =[]
                        dup = tree.get_dup_count()

                        #(this assumes no dups overlap so it's probably overkill 
                        #but should work)
                        for i in dup:   
                            tree.collapse_dups()
                        
                        # Analysing the subtrees. 
                        if tree.istip == True:
                            pass
                        elif tree.label == "D":
                            pass
                        else:
                            conflicts, concordances = comparisons.compare_trees(
                                species_biparts, 
                                species_name_array, 
                                tree, 
                                extra_names[count], 
                                "homolog", 
                                mode, 
                                "some_log_name", 
                                cutoff, 
                                str(file))
                            total_concordances.extend(concordances)

                            filtered_conflicts = comparisons.filter_conflicts(
                            conflicts)
                            total_conflicts.extend(filtered_conflicts)
                                    
                            if ret_perf_concord:
                                if len(filtered_conflicts) != 0:
                                    perfect = False
                            
                            # Storing the subtrees for output. 
                            if outfile_subtrees:
                                tree.get_rid_of_note("CollapsedNotCounted")
                                printout_filename.append(extra_names[count])
                                subtree_printout.append(tree)

            if ret_perf_concord:
                
                # We want the option to filter for only the gene trees that 
                # contain all the species taxa.

                if ret_perf_concord_strict:
                    # Get all the taxa for species and gene trees
                    species_taxa = all_taxa.bipart_proper
                    species_taxa = set(species_taxa)
                    gene_taxa = read_trees.postorder3(gene_root)
                    gene_taxa = gene_taxa.bipart_proper
                    gene_taxa = set(gene_taxa)

                    # Then compare - we want to only include gene trees that 
                    # have all the taxa in the species tree.
                    missing_taxa = species_taxa - gene_taxa
                    if len(missing_taxa) == 0:
                        if perfect:
                            perf_concords.append(file)

                # If filtering not specified it's less complicated
                else:
                    if perfect:
                        perf_concords.append(file)

            # Formatting the subtrees for output.
            if outfile_subtrees:
                outw = open(file+".sub.tre", "w")

                for x in range(0,len(subtree_printout)):
                    subtree_printout[x].get_rid_of_note("CollapsedNotCounted")
                    make_trees.add_loci(subtree_printout[x])
                    
                    if len(printout_filename[x]) != 0:
                        t = "((" + ",".join(printout_filename[x]) + ")-1," + \
                        subtree_printout[x].get_newick_repr() + ");"
                    else:
                        t = subtree_printout[x].get_newick_repr() + ";"
                    
                    # Outputting the subtrees. 
                    outw.write(t+"\n")
                outw.close()                    

        # If we're just looking for the perfectly concordant files, we don't 
        # need to do the rest of the analysis.
        if ret_perf_concord:
            with open(str(outfile_prefix) + "_perfectly_concordant.tsv", 'w') as f:
                for i in perf_concords:
                    f.write(i + "\n")
            sys.exit(0)

        if not args.no_csv:
            # Just writing a more detailed breakdown of the data. 
            outfile=open(str(outfile_prefix) + "_analysis.csv", "w")
            sorted_conflicts=analysis.sort_conflicts(total_conflicts)
            analysis.conflict_stats(sorted_conflicts, species_root, outfile)

        # Making the various output files. 
        concord_out = open(str(outfile_prefix) + "_concord.tre", "w")
        conflict_out = open(str(outfile_prefix) + "_conflict.tre", "w")
        labels_out = open(str(outfile_prefix) + "_labels.tre", "w")
        concord_gene_out = open(str(outfile_prefix) + "_concord_gene_counts.tre", "w")
        conflict_gene_out = open(str(outfile_prefix) + "_conflict_gene_counts.tre", "w")
        if cutoff == 0:
            nodes_analyzable_out = open(str(outfile_prefix) + "_total_analyzed.tre", "w")
        

        # Mapping the concordances and conflicts back onto the species tree. 
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
        
        # Total number of nodes analysed is only output when there is no cutoff
        # because... 
        if cutoff == 0:
            total_rels = []
            total_rels.extend(total_concordances)
            total_rels.extend(total_conflicts)
            make_trees.clear_labels(species_root)
            make_trees.tree_map(species_root, total_rels)
            make_trees.add_zeros(species_root)
            analyzed_tree = species_root.get_newick_repr(showbl=True)
            nodes_analyzable_out.write(analyzed_tree + ";")
            
        # Making the labels tree for referencing the analysis .csv file. 
        make_trees.label_node_ids(species_root)
        labels_tree=species_root.get_newick_repr(showbl=True)
        print "Labels tree: "
        print labels_tree + ";"
        labels_out.write(labels_tree + ";")
    
        # For this, each gene is only counted as concordant or conflicting once 
        # for each node in the species tree. 

        # Concordances:
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
        
        # Conflicts:
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
        
    # Mapping conflicts and concordances onto the gene tree. 
    elif mode == "r":
        
        # Checking for required arguments.
        if not species_tree or not gene_tree:
            print "--species_tree and --gene_tree are required arguments in this \
mode."
            sys.exit()
        
        # Reading in the species tree. 
        file = open(species_tree, "r")
        for line in file:
            tree = line
        
        # Building the species tree and decomposing it into biparts. 
        species_root, species_name_array = make_trees.build(tree)
        species_biparts = read_trees.postorder2(species_root)
        all_taxa = read_trees.postorder3(species_root)
        species_biparts.append(all_taxa)
        
        # Reading in the gene tree and ensuring there is only one. 
        file = open(gene_tree, "r")
        tree = file.readline()
        if file.readline():
            print "There must be only one tree in the gene tree file!"
            sys.exit()
        
        # Building the gene tree from the newick. We make two copies because
        # one of them will be decomposed into subtrees, and end up being 
        # mutated. 
        tree_root, gene_name_array = make_trees.build(tree)
        g, g_array = make_trees.build(tree)
        
        # Mapping the  duplication nodes on the gene tree. 
        make_trees.label_duplications(tree_root)
        make_trees.label_duplications(g)
        
        # Dividing into subtrees. 
        trarray = []
        extra_names = []
        make_trees.subtree_divide(tree_root,trarray, extra_names)
        trarray.append(tree_root)
        extra_names.append([])

        if outfile_subtrees:
            tree_to_print = []
            out_subs = open(gene_tree+".subtree","w")
        
        #parse through the subtrees
        count = -1
        for tree in trarray:
            count += 1

            # Getting the taxa.
            pre_col_taxa = tree.lvsnms()

            # Collapsing duplications to polytomies.
            dup =[]
            dup = tree.get_dup_count()
            for i in dup:   
                tree.collapse_dups()
            
            # If the subtree root is a tip or duplication, it can't be analysed.
            if tree.istip == True:
                pass
            elif tree.label == "D":
                pass

            # All the other subtrees can be analysed. 
            else:
                # Getting the subtree taxa (equivalent to name arrays). 
                subtree_taxa = []
                subtree_taxa = tree.lvsnms()
                subtree_taxa.extend(extra_names[count])
                
                # Bipartitions required to analyse conflict. 
                subtree_biparts = read_trees.postorder2(tree)
                all_taxa = read_trees.postorder3(tree)
                subtree_biparts.append(all_taxa)
                
                # Getting the relationships.
                # Note: cutoff = 0 is always used as "uninformative" is not chosen
                # on that basis.
                rel_list = []
                rel_list = comparisons.comp_biparts(subtree_biparts, species_biparts, subtree_taxa, species_name_array, "some_log_name", 0, "r", "")

                # For the mapping, we want to divide the relationships into
                # conflicts and concordances.
                conflicts = []
                concordances = []

                for rel in rel_list:
                    if rel.relation == 'conflict':
                        conflicts.append(rel)
                    elif rel.relation == 'concordant':
                        concordances.append(rel)
                
                # Locate this subtree in the main tree. 
                mrca = Node()
                g.find_mrca(pre_col_taxa,mrca)
            
                # Labelling the relevant nodes with conflicts and concordances. 
                for x in concordances:
                    mrca.label_node(x.ortholog_bipart,"*",cutoff)
                for x in conflicts:
                    mrca.label_node(x.ortholog_bipart,"X",cutoff)
                
                # Writing out the subtrees if applicable. 
                if outfile_subtrees:
                    out_mrca = copy.deepcopy(tree)
                    out_mrca.clr_label()
                    make_trees.add_loci(out_mrca)
                    if len(extra_names[count]) != 0:
                        out_subs.write("((" + ",".join(extra_names[count]) + ")," + out_mrca.get_newick_repr() + ");\n")
                    else:
                        out_subs.write(out_mrca.get_newick_repr() + ";\n")

        # Adding the loci IDs back to the tree before putting it in the outfile.
        make_trees.add_loci(g)

        # Writing hte labelled tree to the outfile.
        outw = open(str(outfile_prefix) + "_concon.tre", "w")
        outw.write(g.get_newick_repr(showbl=True)+";")
        outw.close()


    # Looking for all gene trees with exact concordance to a given bipart. 
    elif mode == "s": 

        # Checking arguments.
        if not gene_folder or not query_bipart:
            print "--gene_folder and --query_bipart are required arguments in this mode."
            sys.exit()

        # We need this for later.
        concordant_files=[]

        # Putting the query into a form comp_biparts likes.
        query_entered = query_bipart.split(",")

        # Opening and reading the provided folder of homologs.
        if gene_folder[-1] == "/":
            homologs_folder = gene_folder
        else:
            homologs_folder = gene_folder + "/"
        
        file_list = os.listdir(homologs_folder)
        
        # Initialising some things. 
        count = 0
        subtree_printout = []
        printout_filename = []
        final_print_matches = []
        final_print_names = []
        if query_remove:
            removable = query_remove.split(",")

        if outfile_subtrees:
            outw = open(outfile_prefix+".search.subtree","w")

        printout_filename = []
        subtree_printout = []
        names_on_otherside = []
        
        # All files are searched.
        for filename in file_list:
            
            # Housekeeping. 
            count += 1
            total_concordances = []
            sys.stderr.write("processing " + str(count) + " of " + str(len(file_list)) + '\r')

            # Building a tree and making subtrees.
            file_location = str(homologs_folder) + str(filename)
            file = open(file_location, "r")
            tree = file.readline()
            if file.readline():
                print "WARNING!!!! This is only analyzing the first tree in your file"
            
            gene_root, gene_name_array = make_trees.build(tree)
            
            # Labelling duplications. 
            make_trees.label_duplications(gene_root)
            
            # Decomposing the gene tree into subtrees.
            trarray = []
            extra_names = []
            
            make_trees.subtree_divide(gene_root, trarray, extra_names)
            trarray.append(gene_root)
            
            root_names = []
            list_set = set()
            for x in extra_names:
                for y in x:
                    list_set.add(y)
            
            all_names = gene_root.lvsnms()
            root_names = (set(all_names) - list_set)
            extra_names.append(list(root_names))
            
            # Setup.
            # All of the matching biparts will go into ret_bipart.
            ret_bipart = []
            count2 = -1

            for tree in trarray:
                # More setup. 
                count2 += 1
                dup =[]
                dup = tree.get_dup_count()
                
                # Collapsing duplications into polytomies in each subtree.
                for i in dup:   
                    tree.collapse_dups()

                # Neither tips nor duplication nodes are informative.
                if tree.istip == True:
                    pass
                elif tree.label == "D":
                    pass

                # Searching for patterns concordant with the query.
                else:
                    search_list = []
                    
                    # Check if the removable taxa are in this subtree.
                    if query_remove:
                        subtree_taxa = tree.lvsnms()
                        
                        # There cannot be any intersection between the outgroup 
                        # and the query. 
                        check =  any(item in query_entered for item in extra_names[count2])
                        if check:
                            pass
                        
                        # Getting the search list by excluding the removable 
                        # items in the query that are not in this subtree.
                        else:
                            not_in_removable = [] # (but *is* in subtree)
                            not_in_subtree = [] # (but *is* in removable)
                            not_in_removable, not_in_subtree = comparisons.unique_array(removable,subtree_taxa)
                            search_list = [x for x in query_entered if x not in not_in_subtree]
                    else:
                    
                        search_list = query_entered

                    cur_len = len(ret_bipart)

                    # We're only interested if the search query not a tip. NB 
                    # find_bipart expands ret_bipart if it finds a match with 
                    # the search_list in the tree.
                    if 1 < len(search_list):
                        tree.find_bipart(search_list, cutoff, ret_bipart)
                    
                    # If the length of biparts has expanded and we are making a
                    # subtree outfile, get a record of the subtree needed to 
                    # print to the outfile. 
                    if cur_len != len(ret_bipart):
                        if outfile_subtrees:
                            tree.get_rid_of_note("CollapsedNotCounted")
                            printout_filename.append(filename)
                            out_mrca = copy.deepcopy(tree)
                            make_trees.add_loci(out_mrca)
                            subtree_printout.append(out_mrca)
                            names_on_otherside.append(extra_names[count2])

            # Returning the file if we found any matches. 
            if len(ret_bipart) != 0:
                final_print_names.append(filename)
                final_print_matches.append(str(len(ret_bipart)))
                    
        if outfile_subtrees:
            for x in range(0,len(subtree_printout)):
                tree.get_rid_of_note("CollapsedNotCounted")
                make_trees.add_loci(tree)
                t = printout_filename[x] + ": " + subtree_printout[x].get_newick_repr(showbl = True)
                outw.write(t+"\n")
            outw.close()

        # Outputting the summary.
        total_hits = 0
        outw = open(outfile_prefix + "_search.csv", "w")
        for x in range(0, len(final_print_names)):
            outw.write(str(final_print_names[x]) + "," + str(final_print_matches[x])+"\n")
            total_hits += int(final_print_matches[x])
        outw.close()
        print "Total gene families: " + str(len(final_print_names)) + " Total matches: " + str(total_hits)

    # Summarising the output of reverse mode.
    # Outputs totals of each node label.
    # If query_bipart and _species_tree are provided, ...?
    elif mode == "summarize":
        # Check for input.
        if not annotated_tree:
            print "--annotated_tree is required in this mode."
            sys.exit(0)

        # Reading in the annotated tree and bulding it. 
        file = open(annotated_tree, "r")
        tree = file.readline()
        gene_root, gene_name_array = make_trees.build(tree)
        
        # Simple summarising: totals of each node label are printed.
        counts = []
        possible = ["Duplications","D","Concordances","*","Conflicts","X","Uninformative","U"]
        gene_root.count_label(counts)
        for i in range(0,len(possible),2):
            print possible[i] + ": " + str(counts.count(possible[i+1]))
        print "Total tips: " + str(len(gene_name_array))
        
        # What is this....
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
                            outfile_sum.write(str(len(tree.children[0].lvsnms())) \
                                    + "\t" + str(child1_counts.count(possible[1])) \
                                    + "\t" + str(child1_counts.count(possible[3])) \
                                    + "\t" + str(child1_counts.count(possible[5])) \
                                    + "\t" + str(child1_counts.count(possible[7])) \
                                    + "\t" + str(len(tree.children[1].lvsnms())) \
                                    + "\t" + str(child2_counts.count(possible[1])) \
                                    + "\t" + str(child2_counts.count(possible[3])) \
                                    + "\t" + str(child2_counts.count(possible[5])) \
                                    + "\t" + str(child2_counts.count(possible[7])) \
                                    + "\n") 
                    
        else:
            print "For more detailed analysis give species tree and a bipartition of interest"          
            
    #Step1: Get concordance for the orthologs
    #Step2: Get the RT var
    
    elif mode == "sortadate":
		print "Running SortaDate"
		#outgroups = "Struthio_camelus,Tinamou_guttatus"  #Example outgroups to play with later
		if not species_tree or not gene_folder:
			print "--species_tree and --gene_folder are required arguments in \
this mode."
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
		total_concordances = []

        # Making sure the folder name is correct.
		if gene_folder[-1] == '/':
			homologs_folder = gene_folder
		else:
			homologs_folder = gene_folder + '/'

        # Because we do all the work of counting conflicts inside this for loop,
        # we use less memory as we only handle one file at once.
		file_list = os.listdir(homologs_folder)
		len_file_list = str(len(file_list))
		file_no = 0
		sortaout = open(str(outfile_prefix) + "_sorta_date.csv", "w")
		sortaout.write("Filename,Concordances,Root-Tip-Variance,TreeLength,TotalTips\n")

        
        #run through all the files
		for file in file_list:

			file_no += 1
			sys.stderr.write("Processing file " + str(file_no) +
                             " of " + len_file_list + ".\r")

            # Building the gene tree.
			file_location = str(homologs_folder) + str(file)
			gene_file = open(file_location, "r")

            # Reading the gene_file. Can be one or many lines?
			for line in gene_file:

				tree = line
				gene_root, gene_name_array = make_trees.build(tree)
				
				#removable = outgroups.split(",")
				
				#remove the outgroups
				#gene_root = sortadate.remove_ogs(gene_root,removable)
				#gene_root.parent = None

                # This will act on trees with no duplications (orthologs).
				conflicts, concordances = comparisons.compare_trees(species_biparts,
					species_name_array,gene_root,[],
					"ortholog","n","some_log_name",cutoff,str(file))
				
				#Get the RT-Var
				RT_var = sortadate.Get_var(gene_root,gene_name_array)
				
				#Get the treelength
				t_len = sortadate.get_len(gene_root)
				
				sortaout.write(str(file)+","+str(len(concordances))+","+str(RT_var)+","+str(t_len)+","+str(len(gene_name_array))+"\n")				
    else:
        print "Mode not found, for options type: python Conflict_detector.py -h"   
 
