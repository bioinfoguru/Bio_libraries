#!/usr/bin/env python
# coding: utf-8

# ### Motif analysis
# A set of sequences of equal length can be used to create a motif. The `create()` function in the motifs class takes a list of sequences as an argument to create the motif object.

# In[2]:


import Bio
from Bio import motifs
from Bio.Seq import Seq


# In[3]:


instances = [Seq("TACAA"), Seq("TACGC"), Seq("TACAC"), Seq("TACCC"), Seq("AACCC"), Seq("AATGC"), Seq("AATGC")]


# In[4]:


m1 = motifs.create(instances)


# In[5]:


print(m1)


# In[6]:


len(m1)


# The `counts()` function returns the frequency of each alphabet at all the positions in the motif. Frequency of a particular alphabet at different positions within the motif can also be accessed. The `concensus()` function returns the concensus sequence for the motif.

# In[7]:


print(m1.counts)


# In[8]:


m1.counts["A"]


# In[9]:


#get number of "A" at first position
m1.counts["A",0]


# In[10]:


m1.consensus


# To search a particular motif in a sequence, `instances.search()` function can be used.

# In[13]:


sequence_set = [Seq("AACCGGTT"),Seq("AACCCGTT"),Seq("CATTACAA")]
motif_p1 = motifs.create(sequence_set)
motif_p1.instances


# In[74]:


test_seq=Seq("TACACTGCATTACAACCCAAGCATTA")


# In[75]:


for pos, seq in motif_p1.instances.search(test_seq):
    print("%i %s" % (pos, seq))


# ### Motif for protein sequences

# By default, the `create()` function considers the input sequence as a DNA sequence. For creating a motif for protein sequences, the keyword argument `alphabet` need to be specified the all the amino acids in single letter code. 

# In[17]:


instances = [Seq("LXXLL"),Seq("LXXLL")]
motif_prot = motifs.create(instances,alphabet="ACDEFGHIKLMNPQRSTVWXY")


# In[18]:


motif_prot.counts


# ### Motif from MSA
# 

# In[25]:


from Bio import AlignIO


# In[29]:


#Read alignment for proteins with DEAD motif
DEAD_align = AlignIO.read("DEAD2.aln","clustal")
print(DEAD_align)


# In[28]:


print(DEAD_align[1:5,442:448])


# In[46]:


instances_DEAD = [s1.seq for s1 in DEAD_align[:,442:448]]
motif_DEAD = motifs.create(instances_DEAD, alphabet='ACDEFGHIKLMNPQRSTVWY')
#print(motif_DEAD.counts)
df_DEAD_motif = pd.DataFrame.from_dict(motif_DEAD.counts)


# ## Logomaker
# The [logomaker package](https://logomaker.readthedocs.io/en/latest/) offers a rich set of functionality to work with sequences/motifs to create sequence logos. It can be installed via `pip install logomaker`. This library uses pandas and matplotlib to generate sequence logos. The `savefig()` function of the `plt` object can be used to save an image of the logo. The resolution of the reulting image can be adjusted using the `dpi` keyword argument. To draw sequence logo the `Logo()` function can be used which take the sequence motif in the form of pandas dataframe as an argument.

# In[2]:


import logomaker
import pandas as pd
import matplotlib.pyplot as plt


# In[35]:


ss_logo = logomaker.Logo(df_DEAD_motif)
plt.savefig("fig1.png",dpi=600)


# To normalize the values on the y-axis use `normalize()` function.

# In[36]:


motif_pwm = motif_DEAD.counts.normalize()
df_motif_pwm = pd.DataFrame.from_dict(motif_pwm)
ss_logo_pwm = logomaker.Logo(df_motif_pwm)


# ### Decorating logos
# The logo fonts can be changed using the `font_name` argument. Positions within the logo can be highlighted by adding background color as shown below.

# In[84]:


ss_logo_pwm = logomaker.Logo(df_motif_pwm,font_name='Franklin Gothic Book')
ss_logo_pwm.highlight_position(p=1, color='pink')
ss_logo_pwm.highlight_position(p=2, color='pink')
ss_logo_pwm.highlight_position(p=3, color='pink')
ss_logo_pwm.highlight_position(p=4, color='pink')


# ## Reading motifs
# Motif files such available from Jaspar database can be read directly to create a motif object.

# In[49]:


df_jaspar_motif = pd.DataFrame()
fh = open("MA0007.2.jaspar")
for m in motifs.parse(fh, "jaspar"):
    print(m.counts)
    df_jaspar_motif = pd.DataFrame.from_dict(m.counts)
fh.close()


# In[50]:


ss_logo = logomaker.Logo(df_jaspar_motif)


# In[77]:


ss_logo = logomaker.Logo(df_jaspar_motif,font_name='Franklin Gothic Book')
ss_logo.highlight_position(p=3, color='magenta')
ss_logo.highlight_position(p=7, color='lightgreen')

