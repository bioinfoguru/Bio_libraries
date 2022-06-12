#!/usr/bin/env python
# coding: utf-8

# ## Phylogenetic analysis
# One of the frequently used techniques in molecular evolution in phylogenetic analysis. It involves studing the relationship between different organisms based on the similarity between homologous genes. A phylogenetic tree is a pictorial representaion of this relationship. Here we'll use a set of protein sequences to generate a multiple sequence alignment and then based on the distance (which is inverse of similarity) between different sequences we'll construct a phylogenetic tree. There are multiple libraries available to do phylogenetic analysis. We'll use `ete3` library which can be installed using the command `pip install --upgrade ete3`. In addition, to render tree we also need PyQt which can be installed via `pip install PyQt5`.

# In[1]:


from ete3 import Tree


# ### Newick format
# One of the frequently used format for tree representation in Newick format. This is a text based notation of tree structures. In this, nodes are specified as comma separated values within parantheses and ends with a semi-colon. Leaf nodes and internal can have names and branch length can be specified as well. Below are some examples of tree constructed using the `Tree()` function. The `format` attribute of this function need to be changed when constructing tree with flexible newick format. 

# In[2]:


t = Tree( "((a,b),c);" )


# In[3]:


print(t)


# In[26]:


#Tree with internal nodes
t2 = Tree( "((a,b)m,c)n;", format=1 )


# In[27]:


print(t2)


# The `show()` function generates a graphical representation of the tree and `render()` function can used to save an image of the tree. This function has keyword attributes to specify resolution, width and height of the image. To show the tree within notebook we call the `render()` function with `%%inline` argument.

# In[28]:


t.render("tree1.png", dpi=300)
t.render("%%inline")


# To get the names of all the nodes in a tree we can use `traverse()` function to interate through the nodes of the tree and print their names. To get a list of only the leaf nodes, `get_leaves()` can be used.

# In[25]:


print([n.name for n in t2.traverse()])
print([x.name for x in t2.get_leaves()])


# ### Loading free from a file
# We can create a tree object by loading a tree in newick format. In this example we'll use a phylogenetic tree generated after multiple sequence alignment of protein sequences using the MUSCLE program.

# In[30]:


dhfr_tree = Tree('dhfr_20.tree')


# In[31]:


print(dhfr_tree)


# In[32]:


dhfr_tree.render("%%inline")


# ### Decorating trees
# There are function available to add style to the nodes and branches of a tree to highlight specific aspects of the phylogenetic relationship. For example, to highlight a particular clade within the tree, we can add a background color. We can also alter the layout of the tree using the functions available in the `TreeStyle` class. Similarly, nodes can be customized using `NodeStyle`.

# In[33]:


from ete3 import TreeStyle
ts = TreeStyle()
ts.show_leaf_name = True
ts.rotation = 90
dhfr_tree.render("%%inline",tree_style=ts)


# In[34]:


#Circular layout
circular_style = TreeStyle()
circular_style.mode = "c" 
dhfr_tree.render("%%inline", tree_style=circular_style)


# In[35]:


# Semi-circular layout
circular_style.arc_start = -180 
circular_style.arc_span = 180
dhfr_tree.render("%%inline", tree_style=circular_style)


# In[19]:


from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
nst1 = NodeStyle()
nst1["bgcolor"] = "LightSteelBlue"
nst2 = NodeStyle()
nst2["bgcolor"] = "DarkSeaGreen"
n1 = dhfr_tree.get_common_ancestor("sp|Q86XF0|DYR2_HUMAN", "sp|P00374|DYR_HUMAN")
n1.set_style(nst1)
n2 = dhfr_tree.get_common_ancestor("sp|Q04515|DYR10_ECOLX", "sp|P22906|DYR_CANAX")
n2.set_style(nst2)

dhfr_tree.render("%%inline", tree_style=circular_style)


# ### Adding sequence alignment to tree
# For situations where we would like to show sequence alignment (or a subset of alignment) along with the phylogenetic tree, we can use `PhyloTree` class. Here, instead of contructing an object of class `Tree`, the tree in instantiated as an object of the class PhyloTree. This allow to attach a multiple sequence alignment to the phylogenetic tree. The node 

# In[40]:


from ete3 import PhyloTree


# In[42]:


align_tree2 = PhyloTree("( HBAZ_CAPHI:0.36620, HBA_HUMAN:0.07042, HBA_MOUSE:0.07042);")


# In[43]:


align_txt = '''
>HBAZ_CAPHI
MSLTRTERTIILSLWSKISTQADVIGTETLERLFSCYPQAKTYFPHFDLHSGSAQLRAHG
SKVVAAVGDAVKSIDNVTSALSKLSELHAYVLRVDPVNFKFLSHCLLVTLASHFPADFTA
DAHAAWDKFLSIVSGVLTEKYR
>HBA_HUMAN
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP
AVHASLDKFLASVSTVLTSKYR
>HBA_MOUSE
MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHG
KKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTP
AVHASLDKFLASVSTVLTSKYR
'''


# In[44]:


align_tree2.link_to_alignment(alignment=align_txt, alg_format="fasta")


# In[45]:


print(align_tree2)
align_tree2.render("%%inline")
#align_tree2.render(file_name="tree_align.png", dpi=300)

