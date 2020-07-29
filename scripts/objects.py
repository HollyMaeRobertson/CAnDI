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

	def  __init__(self, label = "", unique_id = ""):
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

	def __init__(self, relation = "", species_node = "", ortholog_node = ""):
		self.relation = relation
		self.species_node = species_node
		self.ortholog_node = ortholog_node
		self.species_bipart = []
		self.ortholog_bipart = []

	def add_species_bipart(self, bipart):
		self.species_bipart = bipart
	
	def add_ortholog_bipart(self, bipart):
		self.ortholog_bipart = bipart


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
		self.length = 0.0
		self.time_length = 0.0
		self.parent = None
		self.children = []
		self.data = {}
		self.istip = False
		self.height = 0
		self.locus = ""
		self.unique_id = ""

	def add_child(self, child):
		#make sure that the child is not already in there
		assert child not in self.children
		self.children.append(child)
		child.parent = self
	
	def remove_child(self, child):
		#make sure that the child is in there
		assert child in self.children
		self.children.remove(child)
		child.parent = None
	
	def leaves(self,v=None):
		if v == None:
			v = []
		if len(self.children) == 0:
			v.append(self)
		else:
			for child in self.children:
				child.leaves(v)
		return v

	def leaves_fancy(self):
		return [n for n in self.iternodes() if n.istip ]

	def lvsnms(self):
		return [n.label for n in self.iternodes() if n.istip ]

	def iternodes(self,order="preorder"):
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

	def get_newick_repr(self,showbl=False):
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
