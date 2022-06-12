#!/usr/bin/env python
# coding: utf-8

# ## Introduction to BioPython

# Python offers a variety of functions to work with text data (Strings) that make it easier to work with biological data such as DNA or protein sequences. BioPython library provides a set of classes that dedicated to parsing and analysis of different type of biological data. The functions avaiable in BioPython helps researcher to progammatically process the data. Below we'll see some of the features in Biopython for working with sequence and structure data. To install Biopython library run `pip install biopython`. For more details regarding Biopython installation and tutorials, please refer to the [Biopython wiki](https://biopython.org/wiki/Biopython2).

# To check the version of Biopython, run the following command.

# In[1]:


import Bio
print(Bio.__version__)


# ### Sequence Object
# To work with sequences, we'll need the `Bio.Seq` class which has the required functions for reading and writing sequence data. Once we have imported this class we can create objects having required data. The example below shows constructing a sequence object with a DNA sequence and then using the `complement()` and `translate()` functions to find the sequence of the complementary strand and the translated protein sequence, respectively. 

# In[2]:


from Bio.Seq import Seq
new_seqeunce = Seq('AATTGGAACCTT')
print("Original sequence:", new_seqeunce)
print("The sequence for the complementary strand is: ",new_seqeunce.complement())
print("The translated protein sequence would be: ",new_seqeunce.translate())


# To create a sequence object by reading sequence from a file, we can use the `SeqIO` class. The `parse()` function in this class can read and write sequences in different formats. This function take two arguments - file name and format, and return an iterator having all the sequences. The code below shows reading a file having multiple sequences and printing the sequences using the seq attribute. Note that this would print sequnces without any annotations. The `description` attribute of SeqIO object can used to print the description of a sequence as given in the input file. The sequence file used in the code below is available [here](https://doi.org/10.1371/journal.pcbi.1005975.s002).  

# In[3]:


from Bio import SeqIO
count=0
for all_seqs in SeqIO.parse("kinases.txt", "fasta"):
    print(all_seqs.description)
    print(all_seqs.seq)
    count += 1
    if(count == 3):
        break


# The `write()` function takes three arguments &mdash; 1) a sequence object, 2) filename, and 3) file format. The code below reads a fasta file with multiple sequences and then save the first 10 sequences in a new file. 

# In[12]:


all_seqs = []
for seq_record in SeqIO.parse("kinases.txt", "fasta"):
    all_seqs.append(seq_record)
SeqIO.write(all_seqs[0:10],"test.aln","clustal")


# ### Multiple Sequence Alignment

# The `AlignIO` class has functions to parse alignment files. The `read()` and `write()` functions have a similar syntax to the corresponding functions in the `SeqIO` class. The alignment object stores sequences in 2D array format such that the rows are number of sequences and columns represent alignment length. To extract a sub-set of an alignment, slicing feature can be used. The code below shows reading an alignment file in fasta format followed by selecting a portion of this alignment and save it in a new file in clustal format. The subset is extracted by giving the range for the rows and columns within square brackets. The numbering for both rows and columns starts with zero. In the example below first ten sequences in the alignment are selected since range of rows is `:10` and the colums range is `3:12`. 

# In[16]:


from Bio import AlignIO
align1 = AlignIO.read("kinases.txt", "fasta")
print(align1)


# In[17]:


sub_alignment = align1[:10,3:12] #[row range, column range]
print(sub_alignment)
AlignIO.write(sub_alignment,"msa1.aln","clustal")


# In[ ]:


# %load msa1.aln
CLUSTAL X (1.81) multiple sequence alignment


consensus                           LKVL-GKGA
sp|Q22RR1||agc:agc-sar|Tetrahy      IKEL-GRGN
sp|Q234E6||agc:agc-sar|Tetrahy      IKKL-GFGQ
sp|Q23KG5||agc:agc-sar|Tetrahy      VKKL-GNGQ
sp|Q23DN8||agc:agc-sar|Tetrahy      IKTL-AFGQ
sp|I7MFS4||agc:agc-sar|Tetrahy      IKKL-GVGQ
sp|I7M3B5||agc:agc-sar|Tetrahy      IKKL-GFGQ
sp|I7MD55||agc:agc-sar|Tetrahy      IKKL-GFGQ
sp|Q869J9|pkg-2|agc:pkg|Parame      IKKL-GEGQ
sp|A8N3F0||agc:agc-unique|Copr      IRVL-GKGC



# ## Running BLAST over the internet
# Biopython offers a functionality to programmatically run BLAST on the NCBI servers using the `Bio.Blast` class.

# To run blast online at NCBI servers, `Bio.Blast` can be used which has different function to run Blast and also to parse the output. The `NCBIWWW` library has `qblast()` function that takes three arguments &emdash; 1) blast program (blastp, blastn, etc.), 2) database (any of the databases available at NCBI, and 3) sequence. Once the blast serach is over the output can be saved in a file. This output would be in XML format. You can use `read()` function within the NCBIXML class to parse this output. The code below shows running a blast search using `qblast` against the non-redundant database available at in NCBI.
# The output file saved in the previous step has all the hits identified in the Blast search. These hits follow a hierarchical manner such that each result would have multiple alignments and within each alignment would be multiple high scoring pairs (hsps) i.e. Blast object $\longrightarrow$ Alignment $\longrightarrow$ hsps. For more details on this you may refer to the Blast documentation available at NCBI.
# 

# In[4]:


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

result_ncbi = NCBIWWW.qblast("blastn", "nt", "8332116")
with open("my_blast.xml", "w") as file_handle:
    file_handle.write(result_ncbi.read())


# In[5]:


result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)
count = 0
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print(hsp)
    count += 1
    if(count==5):
        break


# The `hsps` object has several attributes including the Blast statistics such as evalue, score, positives, etc. These can be used to extract hits based on certain conditions. E.g., the code below shows saving hits from the previous Blast search with evalue greater than 1e-105 to a new file. 

# In[21]:


with open("new_file.txt", "w") as file_handle:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if (hsp.expect < 1e-105):
                file_handle.write(str(hsp)+"\n")
            file_handle.write("\n")
print("DONE")


# In[ ]:





# ### BLAST search using sequence file
# To run the Blast search using a sequence file instead of gi number, we first need to create a seqeunce object and then pass it on to the `qblast` function as shown below. To run this code, save the protein sequence below in a new file example1.fasta.

# 
#    MFHPGMTSQPSTSNQMYYDPLYGAEQIVQCNPMDYHQANILCGMQYFNNSHNRYPLLPQMPPQFTNDHPY
#    DFPNVPTISTLDEASSFNGFLIPSQPSSYNNNNISCVFTPTPCTSSQASSQPPPTPTVNPTPIPPNAGAV
#    LTTAMDSCQQISHVLQCYQQGGEDSDFVRKAIESLVKKLKDKRIELDALITAVTSNGKQPTGCVTIQRSL
#    DGRLQVAGRKGVPHVVYARIWRWPKVSKNELVKLVQCQTSSDHPDNICINPYHYERVVSNRITSADQSLH
#    VENSPMKSEYLGDAGVIDSCSDWPNTPPDNNFNGGFAPDQPQLVTPIISDIPIDLNQIYVPTPPQLLDNW
#    CSIIYYELDTPIGETFKVSARDHGKVIVDGGMDPHGENEGRLCLGALSNVHRTEASEKARIHIGRGVELT
#    AHADGNISITSNCKIFVRSGYLDYTHGSEYSSKAHRFTPNESSFTVFDIRWAYMQMLRRSRSSNEAVRAQ
#    AAAVAGYAPMSVMPAIMPDSGVDRMRRDFCTIAISFVKAWGDVYQRKTIKETPCWIEVTLHRPLQILDQL
#    LKNSSQFGSS

# In[14]:


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


# In[16]:


seq_file = open('example1.fasta')
result_handle2 = NCBIWWW.qblast("blastp", "nr", seq_file.read())
seq_file.close()

with open("test_blast.xml", "w") as out_handle:
    out_handle.write(result_handle2.read())

blast_output = open("test_blast.xml")    
blast_record = NCBIXML.read(blast_output)
print(blast_record.alignments[0])


# Let's say we need only the alignment with the mouse sequence, then, to print first 50 characters of each alignment with the mouse sequence along with corresponding statistics, the following code can be used.

# In[17]:


for alignment in blast_record.alignments:
    if "musculus" in alignment.title:
        print(alignment.title)
        for hsp in alignment.hsps:
            print(hsp.query[0:50])
            print(hsp.match[0:50])
            print(hsp.sbjct[0:50])
            print(hsp.positives, hsp.score, hsp.expect)


# To save the Blast output in csv format, we can use the `csv` library as shown below.

# In[19]:


import csv
csv_out = open("blast_out.csv", "w", newline='')
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        new_row = csv.writer(csv_out, delimiter=",")
        new_row.writerow([alignment.title.split("|")[-1], hsp.positives, hsp.expect])
csv_out.close()

