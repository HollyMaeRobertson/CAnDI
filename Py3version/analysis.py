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

        if node not in list(node_dict.keys()):
            node_dict[node] = [conflict]
        else:
            node_dict[node].append(conflict)

    # Now we sort each bipart by type of conflict within its node.
    complete_dict = {}

    for key in list(node_dict.keys()):
        conflict_list = node_dict[key]
        conflict_dict = {}
        conflict_dict["conflict_1"] = [conflict_list[0]]

        for conflict in conflict_list[1:]:
            index = 1
            bp1 = conflict.ortholog_bipart

            while True:
                name = "conflict_" + str(index)

                if name in list(conflict_dict.keys()):
                    bp2 = conflict_dict[name][0].ortholog_bipart
                else:
                    sys.stderr.write("Making list for node_" +
                                     key + " conflict_" + str(index) + "\r")
                    conflict_dict[name] = [conflict]
                    break

                rel = comparisons.bipart_relationship(bp1, bp2)

                if rel == 'concordant':
                    conflict_dict[name].append(conflict)
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
        print("not a list!")


def get_gene_names(con):

    genes = ""
    for i in con:
        genes += ";" + i.gene_name
    return genes[1:]


def conflict_stats(conflicts_dict, tree, outfile):
    """This function should take a dictionary from sort_conflicts and 
    calculate the most common conflict at each node, second-most common, 
    etc.
    """

    # We made this as a dictionary earlier because it was easier to do it
    # that way then, but now we want to put things in a defined order so we
    # need a list.
    stats_dict = {}

    for node in list(conflicts_dict.keys()):
        stats_dict[node] = []

        for name in list(conflicts_dict[node].keys()):
            conflict_list = conflicts_dict[node][name]
            new_list = [name, conflict_list]
            stats_dict[node].append(new_list)

    outfile.write(
        "node_id,species_bipart,ortholog_bipart,alternative_conflicts,number_of_conflicts,percentage,genes\n")

    for node in list(stats_dict.keys()):
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
            output.append(";".join(conflict[1][0].ortholog_bipart))

            # Alternative conflicts should be included where they exist.
            if conflict[1][0].alt_conflict:
                alternatives = []
                alternatives.append(
                    ";".join(sorted(conflict[1][0].alt_conflict)))
                for i in conflict[1]:
                    include = False
                    for j in alternatives:
                        if i.alt_conflict:
                            if ";".join(sorted(i.alt_conflict)) != j:
                                include = True
                    if include:
                        alternatives.append(";".join(sorted(i.alt_conflict)))
                output.append(" : ".join(alternatives))
            else:
                output.append("")
            output.append(str(how_common))
            output.append(str(percent))

            # get the gene names
            gene_names_joined = ""
            gene_names_joined = get_gene_names(conflict[1])
            output.append(gene_names_joined)

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
