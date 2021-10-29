"""This module contains functions for building trees from their newick
representation into a structure that can be traversed easily, using 
Nodes, and functions for editing such trees. 
"""

from objects import Node, Bipart, Rel
import read_trees
import copy


def build(instr):
    """This takes in a tree as a string in newick format, then 
    puts it in a data structure using Nodes that can be easily 
    traversed.
    """
    root = None
    name_array = []
    index = 0
    nextchar = instr[index]
    beginning = True
    keepgoing = True
    current_node = None
    counter = 0

    # The first open bracket in a newick tree format denotes the root node.
    while keepgoing:
        if nextchar == "(" and beginning:
            root = Node()
            current_node = root
            beginning = False

        # In newick format, each new bracket denotes a new node which is
        # the child of the previous node.
        elif nextchar == "(" and not beginning:
            newnode = Node()
            current_node.add_child(newnode)
            current_node = newnode

        # Commas separate taxa in a clade (tips of the tree) so where we
        # see one we should move to a parent node to make the next tip.
        elif nextchar == ',':

            current_node = current_node.parent

        # Closing brackets close clades, so we need to move up a node to
        # open the next one.
        elif nextchar == ")":
            current_node = current_node.parent
            index += 1

            # After closing brackets, there may be labels for that
            # node (e.g. bootstrap values are often found here).
            nextchar = instr[index]

            while True:
                if nextchar == ',' or nextchar == ')' \
                        or nextchar == ':' \
                        or nextchar == ';' \
                        or nextchar == '[':
                    break

                name += nextchar
                index += 1
                nextchar = instr[index]

            #This is to handle floats
            try:
                name = float(name)
                name = int(name)
                name =  str(name)
            except ValueError:
                name = name

            current_node.label = name
            current_node.sup = name

            # We give each node a unique id from the counter, this
            # allows all internal nodes to be easily referenced.
            current_node.unique_id = str(counter)
            counter += 1

            # When we've finished recording the label, we need to go
            # back one to be in the right place for the next step.
            index -= 1

        # In newick format, a semicolon denotes the end of the tree.
        elif nextchar == ';':
            keepgoing = False
            break

        # In newick format, a colon (usually after a label) denotes a
        # branch length, which is information we want to keep.
        elif nextchar == ":":
            index += 1
            nextchar = instr[index]

            while True:
                if nextchar == ',' or nextchar == ')' \
                        or nextchar == ':' \
                        or nextchar == ';' \
                        or nextchar == '[':
                    break

                branch += nextchar
                index += 1
                nextchar = instr[index]

            current_node.length = float(branch)

            # As before, we need to go back a place.
            index -= 1

        # Whitespace means nothing in newick trees.
        elif nextchar == ' ':
            pass

        # Ortholog trees have locus labels preceded by '@'. This locus
        # label shouldn't be part of the taxon label as it refers to a
        # genomic location, not a taxon, but we will still need it.
        elif nextchar == '@':
            while True:
                if nextchar == ',' or nextchar == ')' \
                        or nextchar == ':' \
                        or nextchar == ';' \
                        or nextchar == '[':
                    break

                name += nextchar
                index += 1
                nextchar = instr[index]

            current_node.locus = name
            index -= 1

        # If it's anything else, it's a taxon, so make an external node.
        else:
            newnode = Node()
            current_node.add_child(newnode)
            current_node = newnode
            current_node.istip = True

            while True:
                if nextchar == ',' or nextchar == ')' \
                        or nextchar == ':' \
                        or nextchar == ';' \
                        or nextchar == '[' \
                        or nextchar == '@':
                    break

                name += nextchar
                index += 1
                nextchar = instr[index]

            current_node.label = name
            current_node.sup = name
            name_array.append(name)
            index -= 1

        # Each time, we move on to the next character in the newick.
        if index < len(instr) - 1:
            index += 1
        nextchar = instr[index]

        # Re-initialise these each time.
        name = ""
        branch = ""
        locus = ""

    return root, name_array


#gets the clade that is the first node before the duplication
def get_first_before_dup(root, array):

    #make sure this isn't the overall root
    if root.parent:

        #need the base of the duplication
        if root.parent.label == "D":
            get_first_before_dup(root.parent,array)
        else:
            nd = Node()
            if str(root.parent.children[0].get_newick_repr()) == str(root.get_newick_repr()):
                nd = copy.deepcopy(root.parent.children[1])
                array.append(nd)
            else:
                nd = copy.deepcopy(root.parent.children[0])
                array.append(nd)
            get_first_before_dup(root.parent,array)


#Not a fan of deep copy but this seems to get the job
#done for now
def subtree_divide(root, trarray, extra_names):
    """This takes the root of a gene tree with labelled duplications and a list
    (trarray) which it fills with every node that is a direct child of a 
    duplication node.

    """
    
    # Every child of a duplication gets added to trarray, leading to overlap if 
    # there are repeated duplications.
    if root.label == "D" and root.istip == False:
        Node1 = Node()
        Node2 = Node()
        
        array = []
        name_array = []
        
        # Get the outgroup taxa.
        get_first_before_dup(root,array)
        
        if len(array) != 0:
            for i in array:
                name_array.extend(i.lvsnms_uniq())
        
        # Adding to trarray and extra_names. 
        Node1 = copy.deepcopy(root.children[0])
        Node2 = copy.deepcopy(root.children[1])
        trarray.append(Node1)
        trarray.append(Node2)
        extra_names.append(name_array)
        extra_names.append(name_array)
    
    # Recursion to repeat these steps on every node in the tree. 
    for child in root.children:
        subtree_divide(child, trarray, extra_names)

#Repeat of something not a fan of
def subtree_divide_at_base_of_dup(root, trarray, extra_names):
    
    if root.label == "D" and root.istip == False:
        
        Node1 = Node()
        #need to grab just the clade and make it a polytomy, the tip won't be
        #analyzed as part of a subtree otherwise
        array = []
        name_array = []
        #get the outgroup taxa
        get_first_before_dup(root,array)
        if len(array) != 0:
            for i in array:
                name_array.extend(i.lvsnms_uniq())
        Node1 = copy.deepcopy(root)
        trarray.append(Node1)
        extra_names.append(name_array)
    
    for child in root.children:
        subtree_divide_at_base_of_dup(child,trarray,extra_names)


def add_loci(root):
    # Postorder traversal of tree from specified root, changing all tips
    # into unique values in the form species@locus.

    for child in root.children:
        if child.istip:
            child.label = child.label + child.locus
        add_loci(child)


def subtrees_function(root, subtrees=None):
    """A function to check whether a node is one of various kinds, allowing
    us to check for duplications and hence make the largest possible 
    subtrees without internal duplication events. Takes the root node of the
    tree to be split into subtrees and returns a list of the subtrees. NB:
    the subtree nodes in the list are still the ones in the original tree, 
    this doesn't make copies of them.
    """

    if subtrees is None:
        subtrees = []

    children = [i for i in root.children]

    side1 = children[0]
    side2 = children[1]

    bipart1 = read_trees.postorder3(side1).bipart_proper
    bipart2 = read_trees.postorder3(side2).bipart_proper

    # Checking if *this* node corresponds to a duplication.
    duplication = "No"

    for i in bipart1:
        if i in bipart2:
            duplication = "Yes"
    for i in bipart2:
        if i in bipart1:
            duplication = "Yes"

    # Checking if either side of the node contains another, downstream
    # subduplication.
    side1_dup = len(bipart1) - len(set(bipart1))
    side2_dup = len(bipart2) - len(set(bipart2))

    # If there are subduplications on both sides we need to go further into
    # the tree.
    if side1_dup != 0 and side2_dup != 0:
        subtrees_function(side1, subtrees)
        subtrees_function(side2, subtrees)

    # If only one of the sides contains a subduplication, we make a subtree
    # for the other side and go further into the tree for this side.
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


def label_duplications(root, recursive=True):
    """Labels duplication nodes with 'D'. If recursive is True, traverses 
    the whole tree to do this.
    """
    if root.istip == False:
        children = [i for i in root.children]

        side1 = children[0]
        side2 = children[1]

        if side1.istip:
            bipart1 = [side1.label]
        else:
            bipart1 = read_trees.postorder3(side1).bipart_proper

        if side2.istip:
            bipart2 = [side2.label]
        else:
            bipart2 = read_trees.postorder3(side2).bipart_proper

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


def tree_map(tree_root, bipart_list):
    """This replaces the labels of each node in a species tree with the 
    numbers of entries that apply to that node in bipart_list. 

    WARNING: this function *replaces* the labels in the provided tree, it
    doesn't make a new tree with different ones new ones. Be careful when 
    you call it.
    """

    bipart_dict = {}

    # Count the numbers at each node.
    for bipart in bipart_list:
        key = str(bipart.species_node)
        if key in bipart_dict.keys():
            bipart_dict[key] += 1
        else:
            bipart_dict[key] = 1

    # Add the numbers to the tree. You would not believe how long it took
    # me to get this bit to work for what it is.
    for key in bipart_dict.keys():
        label = str(bipart_dict[key])
        node = read_trees.node_finder(tree_root, key)
        node.label = label


def tree_map2(tree_root, rel_list, label):
    """This replaces the labels of each node in a gene tree with a given
    label, provided they are in rel_list.
    """
    bipart_dict = {}
    if tree_root.parent:
        tree_root = tree_root.parent

    for rel in rel_list:
        key = str(rel.ortholog_node)
        bipart_dict[key] = label

    for key in bipart_dict.keys():
        node = read_trees.node_finder(tree_root, key)
        node.label = label

#Continution the naming scheme
def tree_map3(tree_root, rel_list, outname, multi_info):
    """This replaces the labels of each node in a gene tree with the
    number of times a gene had the relationship
    """
    bipart_dict = {}
    gene_name_dict = {}
    duplicate = {}
    
    #Species node can only informed once per gene tree in this case
    for bipart in rel_list:
        gene_name = str(bipart.gene_name)
        key = str(bipart.species_node)
        
        if key in gene_name_dict.keys():
            if gene_name not in gene_name_dict[key]:
                gene_name_dict[key].append(gene_name)
            #This means the single tree has had more than one con* with the relationship
            else:
                
                #make a key that is both the key and gene name
                tup = (key,gene_name)
                if tup in duplicate.keys():
                    duplicate[tup].append(gene_name)
                #This is the first time it is notices to add it multiple times to account
                #for the first time when it wasn't added
                else:
                    duplicate[tup] = []
                    duplicate[tup].append(gene_name)
                    duplicate[tup].append(gene_name)
        else:
            gene_name_dict[key] = []
            gene_name_dict[key].append(gene_name)

    outw_multi = open(multi_info, "w")
    for key in duplicate.keys():
        outw_multi.write(key[0] + "," + ";".join(duplicate[key])+"\n")
    outw_multi.close()
    
    # Add the numbers to the tree. You would not believe how long it took
    # me to get this bit to work for what it is. (Haha, been there)
    outw = open(outname,"w")
    for key in gene_name_dict.keys():
        outw.write(key + "," + ";".join(gene_name_dict[key])+"\n")
        label = str(len(gene_name_dict[key]))
        node = read_trees.node_finder(tree_root, key)
        node.label = label
    outw.close()

def clear_labels(root):
    """Removes all the labels downstream of the root node specified, except
    taxon names.
    """
    if root.istip == False:
        root.label = ""

    for i in root.children:
        if i.istip == False:
            i.label = ""

        clear_labels(i)


def change_tips_to_species(root):
    """Replaces tip names in the form species@locus with only the species 
    name.
    """
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


def add_zeros(root):
    """Replaces all blank labels in a tree with the label '0'."""
    if root.istip == False:
        if root.label == "":
            root.label = "0"

    for i in root.children:
        if i.istip == False:
            if i.label == "":
                i.label = "0"

        add_zeros(i)


def label_uninformative(root):
    """Replaces all integer labels in a tree with the label 'U'."""
    if root.istip == False:
        if str.isdigit(root.label):
            root.label = "U"

    for i in root.children:
        if i.istip == False:
            if str.isdigit(i.label):
                i.label = "U"

        label_uninformative(i)


def label_node_ids(root):
    """Replaces all labels in a tree with their node id."""
    if root.istip == False:
        root.label = root.unique_id

    for i in root.children:
        if i.istip == False:
            i.label = i.unique_id

        label_node_ids(i)
