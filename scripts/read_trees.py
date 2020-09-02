"""This module contains functions used for 'reading' trees i.e. getting
information from them without changing them. Many of these functions use
recursion to carry out postorder traversals.
"""

from objects import Node, Bipart, Rel
import make_trees


def postorder3(root, bipart=None):
    """This traverses the entire tree downstream of the specified root node
    and puts the names of all the tips into the list 'bipart'.
    """

    # The bipart should always start empty.
    if bipart is None:
        bipart = Bipart(root.label, root.unique_id)

    # We add all the downstream taxa to the bipart.
    for i in root.children:
        if i.istip:
            bipart.add_taxon(i.label)
        postorder3(i, bipart)

    return bipart


def postorder2(root, total_list=None, subtrees=False):
    """This traverses the entire tree downstream of the specified root and 
    makes a list (total_list) of all the possible bipartitions, using 
    postorder3.
    """

    # Total_list should always start off empty
    if total_list is None:
        total_list = []

    # When subtrees is specified as True, we want to make an extra
    # bipartition that encompasses the whole subtree.
    if subtrees:
        bipart = Bipart(root.label, root.unique_id)
        bipart = postorder3(root, bipart)
        total_list.append(bipart)
        subtrees = False

    # Making the rest of the bipartitions.
    for child in root.children:
        if child.children:
            bipart = Bipart(child.label, child.unique_id)

            # Other side of the bipartition (for printing output).
            one_up = child.parent
            side = Bipart()

            for child2 in one_up.children:
                if child2 == child:
                    pass
                else:
                    if child2.istip:
                        side.add_taxon(child2.label)
                    side = postorder3(child2, side)
            bipart.other_side = side.bipart_proper

            # The actual bipart is made and then added to the total
            # list, which we can iterate through
            bipart = postorder3(child, bipart)
            total_list.append(bipart)

        # Recursion so we get all the nodes
        postorder2(child, total_list)

    return total_list


def node_finder(root, node_id):
    """Calls node_finder_recursive to find the node, then returns it -
    required because of the recursive nature of tree traversals.
    """
    node = node_finder_recursive(root, node_id)
    if len(node) == 1:
        node = node[0]
    return node


def node_finder_recursive(root, node_id, node=None, first_time=True):
    """Traverses tree and returns the node with the corresponding id."""

    # Because recursion is annoying, we append our result to an empty list
    # and return the list. This means it keeps its value as it is returned
    # through all the layers.
    if node is None:
        node = []

    if first_time:
        if root.unique_id == node_id:
            node.append(root)

    for i in root.children:
        if not i.istip:
            if i.unique_id == node_id:
                node.append(i)

        node_finder_recursive(i, node_id, node, first_time=False)

    return node


def identify_tricky_nodes(root, subtree_list, tricky_nodes=None, is_root=True):
    """This takes the root node of a tree and a list of all the subtree
    nodes, and returns a list of all the nodes that are not one of the 
    following:
    - a duplication node
    - a subtree root
    - within a subtree.
    """

    # Recursion is a pain so we put things in a list that gets called in every
    # layer.
    if tricky_nodes is None:
        tricky_nodes = []

    # We always check the root first.
    if is_root:
        make_trees.label_duplications(root, recursive=False)

        if root.label != 'D' and root not in subtree_list and root.istip == False:
            tricky_nodes.append(root)
        elif root in subtree_list:
            return tricky_nodes

    # All the other nodes.
    for child in root.children:
        make_trees.label_duplications(child, recursive=False)

        if child.label != 'D' and child not in subtree_list and child.istip == False:
            tricky_nodes.append(child)
        elif child in subtree_list or child.istip == True:
            continue
        identify_tricky_nodes(child, subtree_list, tricky_nodes, is_root=False)

    return tricky_nodes
