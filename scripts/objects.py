"""Key objects for tree traversal.

This module contains the classes Node, Bipart and Rel. Nodes are objects
that describe the nodes on a phylogenetic tree, both internal and
external, while biparts (bipartitions) are objects which correspond to
nodes and allow easy identification of conflict, as they contain all the
taxa downstream of a node.

"""


class Bipart:
    """An object which corresponds to all of the taxa downstream of
    a particular Node (from Node.py). Biparts have the following
    attributes: 
    - unique_id: A unique identifier which corresponds to the Node 
    this bipart is made from.
    - label: The original label on the Node in the tree as entered 
    into the program. It is assumed that the label corresponds to a 
    bootstrap value	or another indicator of confidence in the node, 
    expressed as a 	percentage.
    - other_side: Sometimes empty, other_side is a list of all taxa
    downstream from the 'sister' node of the node this bipart is 
    from.
    - bipart_proper: The bipart proper is a list of all taxa 
    downstream of this bipart's corresponding node. A node can 
    easily be identified on a species tree from this information 
    alone, if necessary.
    """

    def __init__(self, label="", unique_id=""):
        self.label = label
        self.unique_id = unique_id
        self.other_side = []
        self.bipart_proper = []

    def add_taxon(self, taxon):
        self.bipart_proper.append(taxon)

    def add_other_side(self, taxon):
        self.other_side.append(taxon)


class Rel:
    """A class with three main attributes: two node unique_ids and 
    a string containing the relationship between them."""

    def __init__(self, relation="", species_node="", ortholog_node="", gene_name = ""):
        self.relation = relation
        self.species_node = species_node
        self.ortholog_node = ortholog_node
        self.gene_name = gene_name
        self.species_bipart = []
        self.ortholog_bipart = []
        self.alt_conflict = None

    def add_species_bipart(self, bipart):
        self.species_bipart = bipart

    def add_ortholog_bipart(self, bipart):
        self.ortholog_bipart = bipart

    def add_alt_conflict(self, conflict):
        self.alt_conflict = conflict


class Node:
    """An object which corresponds to a node on a phylogenetic tree. 
    Has the following attributes:
    - label: Whatever labels were present on the input tree in 
    newick format - usually assumed to be bootstrap values. String.
    - length: Node 'length' e.g. in substitutions per base pair. 
    Also from the input tree in newick format, often assumed to 
    correspond to evolutionary time. Float.
    - time_length: A variation on length where node_length has been
    corrected for time. Float
    - parent: The parent node of this node, i.e. one node 'upstream'
    on the tree. Node instance.
    - children: The (two) nodes directly downstream of this node. 
    List of node instances.
    - data: I don't know
    - istip: Whether the node in question is a 'tip' i.e. a taxon
    name, or not. Boolean
    - height: I don't know.
    - locus: For gene trees/ortholog trees, the locus is usually a
    part of the label given after the '@' sign in labels of tips. It
    shows the genomic locus of the ortholog, or some other unique 
    genomic identifier of that gene. String.
    - unique_id: A unique node identifier, as a number, that ensures 
    biparts can always be related back to nodes. String.
    """

    def __init__(self):
        self.label = ""
        self.sup = "" #A support that should not be changed 
        self.length = 0.0
        self.time_length = 0.0
        self.parent = None
        self.children = []
        self.data = {}
        self.istip = False
        self.height = 0
        self.locus = ""
        self.unique_id = ""

    #clears info placed on internal nodes
    def clr_label(self):
    
        if self.istip == False and self.label != "D":
            self.label = self.sup
        for child in self.children:
            child.clr_label()
    
    def get_rid_of_note(self,note):
    
        if self.label == note:
            self.label = "-1"
        for child in self.children:
            child.get_rid_of_note(note)
    
    
    #searches for an exact match of a bipartition in question
    def find_bipart(self,bipart_to_find,cutoff,ret_bip):
        
        #get bipart
        bipart = self.lvsnms()
        if sorted(bipart) == sorted(bipart_to_find) and self.label != "CollapsedNotCounted":
            if self.label:
                if str.isdigit(self.label):
                    if int(self.label) < int(cutoff):
                        pass
                    else:
                        ret_bip.append(bipart)
            else:
                ret_bip.append(bipart)
        
        for child in self.children:
        	child.find_bipart(bipart_to_find,cutoff,ret_bip)
    
    #Labels nodes as with conflict or a U assuming they are below a certain
    #cutoff
    def label_node(self,rel,sign,cutoff):
        
        array = []
        array = self.lvsnms_uniq()

        if sorted(array) == sorted(rel):

            #check if there is a label
        	if self.label:
        	    #make sure the label is a digit
        	    if str.isdigit(self.label):
        	        #Anything below a certain level should be labeled U
        	        if int(self.label) < int(cutoff):
        	            self.label = "U"
        	        else:
        	            self.label = sign
        	else:
        		self.label = sign

        for clade in self.children:
            clade.label_node(rel,sign,cutoff)
       
    
    #Finds and adds to an mrca node structure the mrca of
    #a list of taxa
    def find_mrca(self,taxa_list,mrca):

        array = []
        array = self.lvsnms()
        if sorted(array) == sorted(taxa_list):
        	#print array
        	#mrca.children.append(self.children[0])
        	#mrca.children.append(self.children[1])
        	mrca.add_child(self)
        	#print mrca.get_newick_repr()
        
        #allows for postorder without iternodes because why not
        for clade in self.children:          
            clade.find_mrca(taxa_list,mrca)
    
    
    def add_child(self, child):
        # make sure that the child is not already in there
        assert child not in self.children
        self.children.append(child)
        child.parent = self

    def remove_child(self, child):
        # make sure that the child is in there
        assert child in self.children
        self.children.remove(child)
        child.parent = None

    def leaves(self, v=None):
        if v == None:
            v = []
        if len(self.children) == 0:
            v.append(self)
        else:
            for child in self.children:
                child.leaves(v)
        return v

    def leaves_fancy(self):
        return [n for n in self.iternodes() if n.istip]

    def lvsnms(self):
		return [n.label for n in self.iternodes() if n.istip]
	
    def lvsnms_uniq(self):
        array = [n.label for n in self.iternodes() if n.istip]
        set_list = set(array)
        u_list = list(set_list)
        return u_list
            
    def get_dup_count(self):
        return [n.label for n in self.iternodes() if n.label == "D"]
    
	#Turns an edge into a polytomy
    def collapse_node(self,child,Message):
    	
    	for clade in self.children:
    	   
    	   #Not a fan of this but it seems vigorous so far, need to 
    	   #rewrite at some point
    	   if child.get_newick_repr() == clade.get_newick_repr():
               nd = Node()
    	       nd = child.collapse()
    	       self.remove_child(clade)
    	       #if self.parent != None:
    	       #    self.remove_kink(child)
    	       nd.label = Message
    	       self.add_child(nd)
           else:
               clade.collapse_node(child,Message)

    
    #collapses duplication and puts a note to not use the polytomy in downstream analysis
    def collapse_dups(self):
    	
    	for child in self.children:
    		if child.label == "D":
    			#if a duplication only has one child on a side and the parent is a tip
    			#this allows it to be evaluated
    			if child.parent.children[0].istip or child.parent.children[1].istip:
    				self.collapse_node(child,"CollapsedNotCounted")

    				#Account for a bunch of dups in a row
    				#Also ensures they don't lose support value of bipartition
    				if self.parent:
    				    if self.label:
    					    self.parent.collapse_node(self,self.label)
    				    else:
    					    self.parent.collapse_node(self,"")
    			else:
    				self.collapse_node(child,"CollapsedNotCounted")
    		child.collapse_dups()
    		
    def count_label(self,counts):
    	
    	if self.istip == False:
    		counts.append(self.label)
    	
    	for child in self.children:
    		child.count_label(counts)
    
    #Turns the whole tree into a polytomy  
    def collapse(self):
        array = self.lvsnms_uniq()
        main_t = Node()
        for i in array:
            nd = Node()
            nd.label = i
            nd.istip= True
            main_t.add_child(nd)		   
        return main_t		
	
	
	#Removes anything following the duplication	   
    def lvsnms1(self):
        
        for child in self.children:
        	
        	if child.label == "D":

        		self.remove_child(child)
        		if self.parent != None:
        			self.remove_kink(child)
        		main_tree = Node()	
        	child.lvsnms1()

	
	'''
	Modified from PyPHLAWD, kink is an extra edge
	that is left over from removing a child, this removes
	that extra edge	
	'''
    def remove_kink(self,curroot):
        """
        smooth the kink created by prunning
        to prevent creating orphaned tips
        after prunning twice at the same node
        """
        if self == curroot: #and len(curroot.children) == 2:
            #move the root away to an adjacent none-tip
            if curroot.children[0].istip: #the other child is not tip
                curroot = reroot(curroot,curroot.children[1])
            else: curroot = reroot(curroot,curroot.children[0])
        #---node---< all nodes should have one child only now
        length = self.length + (self.children[0]).length
        par = self.parent
        kink = self
        self = self.children[0]
        #parent--kink---node<
        par.remove_child(kink)
        par.add_child(self)
        self.length = length
    
    def iternodes(self, order="preorder"):
        if order.lower() == "preorder":
            yield self
        for child in self.children:
            for d in child.iternodes(order):
                yield d
        if order.lower() == "postorder":
            yield self
	
    def prune(self):
        p = self.parent
        if p != None:
            p.remove_child(self)
        return p

    def get_newick_repr(self, showbl=False):
        ret = ""
        for i in range(len(self.children)):
            if i == 0:
                ret += "("
            ret += self.children[i].get_newick_repr(showbl)
            if i == len(self.children)-1:
                ret += ")"
            else:
                ret += ","
        if self.label != None:
            ret += str(self.label)
        if showbl == True:
            ret += ":" + str(self.length)
        return ret
	
    def set_height(self):
        if len(self.children) == 0:
            self.height = 0
        else:
            tnode = self
            h = 0
            while len(tnode.children) > 0:
                if tnode.children[1].length < tnode.children[0].length:
                    tnode = tnode.children[1]
                else:
                    tnode = tnode.children[0]
                h += tnode.length
            self.height = h
