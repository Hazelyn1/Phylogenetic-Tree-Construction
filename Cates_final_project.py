#Hazelyn Cates
#Started 11/6/23
#EN.605.651.81.FA23
#This program constructs a phylogenetic tree using Biopython package

#DELETE THE FOLLOWING COMMENT, ADD TO REFERENCES IN REPORT!!!!
#Largely based on this tutorial: https://medium.com/@poudelmohit59/beginners-guide-to-phylogenetic-tree-construction-using-biopython-5accbd8345a2

import Bio
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

#read in FASTA sequences from file ad put them in an array that can be indexed:

#10 FASTA sequences were aligned via the online tool MUSCLE (https://www.ebi.ac.uk/Tools/msa/muscle/)
#The alignment file (in .clw format) is read in using AlignIO package:
with open("msa_muscle_output.clw", "r") as msa:
    aligned_seqs = AlignIO.read(msa, "clustal")

print(aligned_seqs)

#FIX THIS!!! Need to somehow change the names of each sequence b/c it's just the accession numbers rn and I have no idea which is which
#But if I can't, I'll make sure I specify what's what in my report when I do my results section

#From the alignment data, have to create a distance matrix
#To do this, first calculate the distances between the sequences: (NOTE: I'm at part 3 in the tutorial linked above)
calc = DistanceCalculator('identity') #for DNA distance calculations

#Now generate the distance matrix of the aligned sequences using this distance calculator
dist_matrix = calc.get_distance(aligned_seqs)
#print(dist_matrix)

#From here, the NJ and UPGMA algorithms can be used to generate trees and compare the results
#First, the NJ method: https://homolog.us/Biopython/Bio.Phylo.TreeConstruction.DistanceTreeConstructor.html#content
construct = DistanceTreeConstructor()
NJ_tree = construct.nj(dist_matrix)
#Phylo.draw(NJ_tree)

#Repeat for UPGMA method:
UPGMA_tree = construct.upgma(dist_matrix)

#print(UPGMA_tree)
#Phylo.draw(UPGMA_tree)
