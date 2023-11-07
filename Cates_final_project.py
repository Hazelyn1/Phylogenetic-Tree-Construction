#Hazelyn Cates
#Started 11/6/23
#EN.605.651.81.FA23
#This program constructs a phylogenetic tree

#DELETE THE FOLLOWING COMMENT, ADD TO REFERENCES IN REPORT!!!!
#Largely based on this tutorial: https://medium.com/@poudelmohit59/beginners-guide-to-phylogenetic-tree-construction-using-biopython-5accbd8345a2

import Bio
from Bio import SeqIO, AlignIO, Phylo

#read in FASTA sequences from file ad put them in an array that can be indexed:

#10 FASTA sequences were aligned via the online tool MUSCLE (https://www.ebi.ac.uk/Tools/msa/muscle/)
#The alignment file (in .clw format) is read in using AlignIO package:
with open("msa_muscle_output.clw", "r") as msa:
    aligned_seqs = AlignIO.read(msa, "clustal")

print(aligned_seqs)

#FIX THIS!!! Need to somehow change the names of each sequence b/c it's just the accession numbers rn and I have no idea which is which
#But if I can't, I'll make sure I specify what's what in my report when I do my results section

#From the alignment data, have to create a distance matrix
#To do this, first calculate the distances between the sequences:

#NOTE: I'm at part 3 in the tutorial linked above
