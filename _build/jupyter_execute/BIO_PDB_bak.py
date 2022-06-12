#!/usr/bin/env python
# coding: utf-8

# ## Working with PDB files

# ### PDB format
# Biomolecular structures can be saved in different formats such as PDB, mmCIF, XML, etc. The PDB format is the most common format and is supported by all visualization and docking programs. The PDB class of Bio library has a `PDBParser()` function to instantiate a PDB parser object. 

# In[ ]:


# %load 7jtl_atoms.pdb
HEADER    VIRAL PROTEIN                           18-AUG-20   7JTL              
TITLE     STRUCTURE OF SARS-COV-2 ORF8 ACCESSORY PROTEIN                        
COMPND    MOL_ID: 1;                                                            
COMPND   2 MOLECULE: ORF8 PROTEIN;                                              
COMPND   3 CHAIN: A, B;                                                         
COMPND   4 SYNONYM: ORF8, NON-STRUCTURAL PROTEIN 8, NS8;                        
COMPND   5 ENGINEERED: YES                                                      
SOURCE    MOL_ID: 1;                                                            
SOURCE   2 ORGANISM_SCIENTIFIC: SEVERE ACUTE RESPIRATORY SYNDROME CORONAVIRUS   
SOURCE   3 2;                                                                   
SOURCE   4 ORGANISM_COMMON: 2019-NCOV;                                          
SOURCE   5 ORGANISM_TAXID: 2697049;                                             
SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;                                 
SOURCE   7 EXPRESSION_SYSTEM_TAXID: 562                                         
ATOM      1  N   GLN A  18      40.917  33.173 126.601  1.00 73.92           N  
ATOM      2  CA  GLN A  18      42.025  33.162 125.645  1.00 73.89           C  
ATOM      3  C   GLN A  18      41.754  32.204 124.493  1.00 67.88           C  
ATOM      4  O   GLN A  18      42.414  31.177 124.374  1.00 66.48           O  
ATOM      5  CB  GLN A  18      43.339  32.782 126.337  1.00 73.82           C  
ATOM      6  N   GLU A  19      40.764  32.527 123.665  1.00 60.90           N  
ATOM      7  CA  GLU A  19      40.434  31.695 122.517  1.00 50.80           C  
ATOM      8  C   GLU A  19      40.928  32.273 121.201  1.00 50.81           C  
ATOM      9  O   GLU A  19      41.122  31.517 120.240  1.00 47.32           O  
ATOM     10  CB  GLU A  19      38.918  31.467 122.431  1.00 55.27           C  
ATOM     11  N   CYS A  20      41.101  33.590 121.115  1.00 51.03           N  
ATOM     12  CA  CYS A  20      41.551  34.244 119.893  1.00 41.77           C  
ATOM     13  C   CYS A  20      42.699  35.209 120.189  1.00 40.41           C  
ATOM     14  O   CYS A  20      42.709  35.887 121.221  1.00 46.45           O  
ATOM     15  CB  CYS A  20      40.393  34.986 119.225  1.00 46.48           C  
ATOM     16  SG  CYS A  20      40.867  35.722 117.656  1.00 56.21           S  
ATOM     17  N   SER A  21      43.669  35.254 119.283  1.00 38.94           N  
ATOM     18  CA  SER A  21      44.851  36.104 119.410  1.00 35.69           C  
ATOM     19  C   SER A  21      44.916  36.902 118.118  1.00 35.79           C  
ATOM     20  O   SER A  21      45.172  36.333 117.051  1.00 32.57           O  
ATOM     21  CB  SER A  21      46.118  35.269 119.634  1.00 33.65           C  
ATOM     22  OG  SER A  21      47.328  36.023 119.512  1.00 33.95           O  
ATOM     23  N   LEU A  22      44.643  38.202 118.205  1.00 36.94           N  
ATOM     24  CA  LEU A  22      44.631  39.082 117.041  1.00 35.49           C  
ATOM     25  C   LEU A  22      45.993  39.763 116.942  1.00 34.55           C  
ATOM     26  O   LEU A  22      46.412  40.451 117.875  1.00 37.54           O  
ATOM     27  CB  LEU A  22      43.520  40.122 117.155  1.00 38.67           C  
ATOM     28  CG  LEU A  22      43.475  41.147 116.024  1.00 41.31           C  
ATOM     29  CD1 LEU A  22      43.308  40.438 114.693  1.00 37.10           C  
ATOM     30  CD2 LEU A  22      42.386  42.162 116.241  1.00 46.10           C  
ATOM     31  N   GLN A  23      46.676  39.584 115.820  1.00 32.31           N  
ATOM     32  CA  GLN A  23      48.026  40.101 115.679  1.00 35.96           C  
ATOM     33  C   GLN A  23      48.129  40.869 114.379  1.00 36.92           C  
ATOM     34  O   GLN A  23      47.308  40.706 113.477  1.00 34.85           O  
ATOM     35  CB  GLN A  23      49.080  38.978 115.752  1.00 33.22           C  
ATOM     36  CG  GLN A  23      49.081  38.307 117.122  1.00 36.84           C  
ATOM     37  CD  GLN A  23      50.113  37.220 117.258  1.00 37.90           C  
ATOM     38  OE1 GLN A  23      51.178  37.277 116.645  1.00 33.05           O  
ATOM     39  NE2 GLN A  23      49.808  36.219 118.072  1.00 35.01           N  
ATOM     40  N   SER A  24      49.141  41.729 114.310  1.00 33.69           N  
ATOM     41  CA  SER A  24      49.355  42.636 113.187  1.00 36.00           C  
ATOM     42  C   SER A  24      50.733  42.408 112.590  1.00 38.21           C  
ATOM     43  O   SER A  24      51.714  42.271 113.320  1.00 33.70           O  
ATOM     44  CB  SER A  24      49.242  44.119 113.611  1.00 38.78           C  
ATOM     45  OG  SER A  24      47.976  44.409 114.168  1.00 44.24           O  


# In[1]:


from Bio.PDB import *


# In[2]:


parser = PDBParser()
structure = parser.get_structure("Protease", "7jtl.pdb")


# In[3]:


print(structure.header["name"])
print(structure.header["resolution"])
print(structure.header["keywords"])


# ## Structure object
# Structure $\longrightarrow$ Model $\longrightarrow$ Chain $\longrightarrow$ Residue $\longrightarrow$ Atom

# In[4]:


#for model in structure:
#    print(model)
for model in structure:
    for chain in model:
        print(chain)


# In[8]:


ctr=0
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if(ctr==5):
                    break
                print(atom)
                ctr+=1


# In[27]:


ctr=0
prot_seq=""
for residue in model.get_residues():
    if(ctr==5):
        break
    print(residue.resname, residue.id[1])
    prot_seq+=residue.resname
    ctr+=1
print(prot_seq)


# To convert the sequence in three letter code to single letter code, we can use `SeqUtils` package available in biopython. The `seq1()` function returns single letter sequence for a three letter sequence. Similarly, to convert single letter sequence to three letter we have `seq3()` function. 

# In[28]:


from Bio.SeqUtils import seq1
seq1(prot_seq)


# In case the structure has some modified residues then custom mapping of three letter code to single letter code can be done using the `custom_map` keyword argument to `seq1()` function. E.g., let say a protein has a phosphorylated serine residue for which the three letter code is SEP and there is no standard single letter code for this. We can add custom mapping for SEP to S as follows.

# In[29]:


prot_seq+="SEP"
print(seq1(prot_seq)) #without custom map
print(seq1(prot_seq,custom_map={"SEP": "S"}))


# ### Shortcuts for accessing structure elements

# In[38]:


atom_list = Selection.unfold_entities(structure, "A")
all_coords = [a.coord for a in atom_list]
print(all_coords[0:3])


# #### Nested lists
# 

# In[28]:


list1 = [[[1,2],[4,5]],[[6,7],[8,9]]]


# In[29]:


list1[0][0][1]


# In[ ]:


#atom = structure[0]["A"][100]["CA"]
#atom


# ### Residue object

# In[39]:


res_list = Selection.unfold_entities(structure, "R")
print(res_list[3].get_resname())
print(res_list[3].get_list())
print(res_list[3].get_id()[1])


# ### Atom object

# In[40]:


atoms_resi3 = res_list[3].get_list()
print(atoms_resi3)
print(atoms_resi3[1].get_name())
print(atoms_resi3[1].get_coord())


# ### Distance Calculations

# In[42]:


ca_3 = res_list[3].get_list()[1]
ca_15 = res_list[15].get_list()[1]


# In[43]:


distance = ca_3 - ca_15
print(distance)


# Calulating Ca-Ca distance matrix

# In[70]:


atom_list = res_list = Selection.unfold_entities(structure, "A")
ca_ca_dist = []

print(len(atom_list))
ca_atom_list = [x for x in atom_list if "CA" in x.name]
print(len(ca_atom_list))

for x in range(len(ca_atom_list)):
    for y in range(x + 1,len(ca_atom_list)):
        ca_ca_dist.append(ca_atom_list[x]-ca_atom_list[y])
print(len(ca_ca_dist))


# In[79]:


from scipy.spatial.distance import squareform
dist_mat = squareform(ca_ca_dist)


# In[84]:


print(dist_mat.shape)
dist_mat


# In[114]:


import matplotlib.pyplot as plt

plt.imshow(dist_mat, cmap="RdBu_r")
plt.colorbar(label="Distance ($\AA$)")
plt.show()


# In[126]:


cont_map = dist_mat
cont_map[cont_map>7]=0


# In[132]:


plt.imshow(cont_map, cmap="YlGn")
plt.colorbar(label="Distance ($\AA$)")
plt.show()


# In[14]:


A_atoms = Selection.unfold_entities(structure[0]["A"], "A") #entity list, target level (A,R,C,M,S)
B_atoms = Selection.unfold_entities(structure[0]["B"], "A")


# In[ ]:





# In[ ]:





# In[ ]:





# In[15]:


ns = NeighborSearch(A_atoms)

contacts_B = []
for atom in B_atoms:
    close_atoms = ns.search(atom.coord, 3)
    if len(close_atoms)>0:
        for a in close_atoms:
            contacts_B.append(a.get_parent().get_id()[1])


# In[16]:


set(contacts_B)


# In[7]:


ns = NeighborSearch(B_atoms)

contacts_A = []
for atom in A_atoms:
    close_atoms = ns.search(atom.coord, 3)
    if len(close_atoms)>0:
        for a in close_atoms:
            contacts_A.append(a.get_parent().get_id()[1])


# In[8]:


set(contacts_A)


# In[ ]:





# ### RMSD

# In[43]:


import numpy as np


# In[44]:


from Bio.PDB.QCPSuperimposer import *
sp1 = QCPSuperimposer()


# In[45]:


print(sp1)


# In[55]:


sp1.set(np.array(all_coords),np.array(all_coords))


# In[53]:


sp1.rms


# In[1]:


import nglview


# In[6]:


view = nglview.show_structure_file("7jtl.pdb")
view.add_representation('cartoon', selection='protein', color='blue')
view


# In[5]:


print(type(view))


# ## Binana

# In[60]:


import sys
sys.path.append("/Users/manishdatt/Desktop/Workshop_NGH/Python/workshop_material/binana-2.1/python/")
import binana


# In[62]:


ligand, receptor = binana.load_ligand_receptor.from_files("/Users/manishdatt/Desktop/Workshop_NGH/Python/workshop_material/binana-2.1/python/example/ligand.pdbqt",                                                           "/Users/manishdatt/Desktop/Workshop_NGH/Python/workshop_material/binana-2.1/python/example/receptor.pdbqt")


# In[63]:


hbond_inf = binana.interactions.get_hydrogen_bonds(ligand, receptor)


# In[64]:


for hbond_label in hbond_inf["labels"]:
    print(hbond_label)


# In[68]:


all_inf = binana.interactions.get_all_interactions(ligand, receptor) 


# In[66]:


all_inf.keys()


# In[67]:


all_inf["hydrogen_bonds"]


# In[69]:


all_data = binana.output.dictionary.collect_all(all_inf)


# In[73]:


print(len(all_data["hydrogenBonds"]))
for x in all_data["hydrogenBonds"]:
    print(x)


# In[4]:





# In[3]:


import ipywidgets 
ipywidgets.Text("hello") 

