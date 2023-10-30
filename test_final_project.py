#Hazelyn Cates
#10/26/23
#Testing stuff for final project on phylogenetic tree construction

from Bio import Phylo, AlignIO, SeqIO

"""fasta_sequences = SeqIO.parse(open(input_file),'fasta')
with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        new_sequence = some_function(sequence)
        write_fasta(out_file)"""

#This is my code, adapted from the above code
"""sequences = SeqIO.parse(open("sequences.txt"), 'fasta')
with open(sequences) as seq_output:
    for f in sequences:
        name, sequence = f.id, str(f.seq)
        write_fasta(seq_output)
"""

#Rad this works (https://biopython.org/wiki/SeqIO)
for line in SeqIO.parse("aligned_seqs.fas", "fasta"):
    print(line)
#Problem is, this doesn't save the file contents to anything, it just prints them out

#This don't work...need to work through this tutorial:
#https://medium.com/@poudelmohit59/beginners-guide-to-phylogenetic-tree-construction-using-biopython-5accbd8345a2
with open("aligned_seqs.fas", "r") as f:
    alignment = AlignIO.read(f, "clustal")

print(alignment)



#So that's my alignment file, so I don't need to do the alignment I just need to build the tree now


#But I'm going to start with aligning like 4-5 FASTA sequences of the TNF gene from different organisms
#So I can feed that alignment information in and make a phylogenetic tree

#Have to perform sequence alignment first, which I can do using AlignIO in biopython OR
#Using clustal omega and importing the results file (see "medium" link above)
#Second, I calculate a distance matrix of the sequences using DistanceCalculator also in biopython
#Then I can use the DistanceTreeConstructor function to build the phylogenetic tree specifying the method (https://biopython.org/wiki/Phylo)

#Seems simple enough, but I can already tell it's going to be a pain in my ass lol 


