#Hazelyn Cates
#Started 11/6/23, finished 11/16/23
#EN.605.651.81.FA23
#This program constructs a phylogenetic tree using Biopython

import Bio
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def four_point_condition(dist_matrix): #checks if the calculated distance matrix is additive
    from itertools import permutations
    import math

    flag = 0 #flag to keep track of if matrix is additive. 0 = not additive
    #First, calculate how many combinations of 4 there are:
    print("Enter number of sequences:")
    n = int(input())
    k = 4

    num_combos = math.comb(n, k) 
    #print(num_combos)

    names = dist_matrix.names #accession numbers

    #print(names.index("NM_001081819.2"))
    #Now find all the combos of 4:
    perms = permutations(names, 4)
    combos = [] #array to hold the 210 permutations

    for i in list(perms): #add all permutations to a list
        #print(i)
        combos.append(i)
    #rint(combos)

    for i in range(0, num_combos): #number of possible permutations of the 10 sequences (210 total)
        a = names.index(combos[i][0]) #i = index of permutation sequence with the "combos" list (it's a list of lists)
        b = names.index(combos[i][1])
        c = names.index(combos[i][2])
        d = names.index(combos[i][3])

        #Now I need to apply the following formula
        #max{𝑑(𝑎,𝑐)+𝑑(𝑏,𝑑),𝑑(𝑎,𝑑)+𝑑(𝑏,𝑐)}≥𝑑(𝑎,𝑏)+𝑑(𝑑,𝑐) --> see M9 assignment & M9 notes for reference
        #Where a = 1, b = 2, c = 3, d = 4
        #So when accessing these combination, need to access a positin WITHIN the list (is that possible?? idk)
        #So like in the distance matrix, I need to access the correct values given the accession numbers in the current combo of 4

        val_1 = max((round(dist_matrix[a][c], 1) + round(dist_matrix[b][d], 1), round(dist_matrix[a][d],1) + round(dist_matrix[b][c], 1)))
        #print(val_1)
        val_2 = round(dist_matrix[a][b], 1) + round(dist_matrix[d][c], 1)
        #print(val_2)
        if val_1 >= val_2: #So far, the matrix is additive
            flag = 1
            continue

        else: #matrix is not additive
            flag = 0
            break

    if flag == 1:
        print("Matrix is additive\n")
    else:
        print("Matrix is not additive\n")

#FASTA sequences were aligned via the online tool MUSCLE (https://www.ebi.ac.uk/Tools/msa/muscle/)
#The alignment file (in .clw format) is read in using AlignIO package:
file_name = input("Enter .clw file name (please include file extension in name)\n") #ex.) msa_muscle_output.clw

with open(file_name, "r") as msa:
    aligned_seqs = AlignIO.read(msa, "clustal")

#From the alignment data, have to create a distance matrix
#To do this, first calculate the distances between the sequences: (NOTE: I'm at part 3 in the tutorial linked above)
calc = DistanceCalculator('identity') #for DNA distance calculations

#Now generate the distance matrix of the aligned sequences using this distance calculator
dist_matrix = calc.get_distance(aligned_seqs)

#The following function does that, but needs to be COMPLETED!
four_point_condition(dist_matrix)
#print(dist_matrix)


#From here, the NJ and UPGMA algorithms can be used to generate trees and compare the results
#First, the NJ method: https://homolog.us/Biopython/Bio.Phylo.TreeConstruction.DistanceTreeConstructor.html#content
construct = DistanceTreeConstructor()
NJ_tree = construct.nj(dist_matrix)
#Phylo.draw(NJ_tree)

#Repeat for UPGMA method:
UPGMA_tree = construct.upgma(dist_matrix)

print("Print NJ tree - N\nPrint UPGMA tree - U\nPrint both trees - B")
tree_input = str(input())
tree_input = tree_input.upper()

if tree_input == 'N': #Draw NJ tree
    Phylo.draw(NJ_tree)
    print("Branch lengths of terminal nodes (leaves) of NJ tree:")
    print(NJ_tree.get_terminals())
    print("Branch lengths of non-terminal (internal) nodes of NJ tree:")
    print(NJ_tree.get_nonterminals())
    print("Print UPGMA tree? Y or N")
    ans = str(input()).upper()
    if ans == 'Y':
        Phylo.draw(UPGMA_tree)
        print("Branch lengths of terminal nodes (leaves) of UPGMA tree:")
        print(UPGMA_tree.get_terminals())
        print("Branch lengths of non-terminal (internal) nodes of UPGMA tree:")
        print(UPGMA_tree.get_nonterminals())

elif tree_input == 'U': #Draw UPGMA tree
    Phylo.draw(UPGMA_tree)
    print("Branch lengths of terminal nodes (leaves) of UPGMA tree:")
    print(UPGMA_tree.get_terminals())
    print("Branch lengths of non-terminal (internal) nodes of UPGMA tree:")
    print(UPGMA_tree.get_nonterminals())
    print("Print NJ tree? Y or N")
    ans = str(input()).upper()
    if ans == 'Y':
        Phylo.draw(NJ_tree)
        print("Branch lengths of terminal nodes (leaves) of NJ tree:")
        print(NJ_tree.get_terminals())
        print("Branch lengths of non-terminal (internal) nodes of NJ tree:")
        print(NJ_tree.get_nonterminals())

else: #Draw both
    Phylo.draw(NJ_tree)
    print("Branch lengths of terminal nodes (leaves) of NJ tree:")
    print(NJ_tree.get_terminals())
    print("Branch lengths of non-terminal (internal) nodes of NJ tree:")
    print(NJ_tree.get_nonterminals())
    Phylo.draw(UPGMA_tree)
    print("Branch lengths of terminal nodes (leaves) of UPGMA tree:")
    print(UPGMA_tree.get_terminals())
    print("Branch lengths of non-terminal (internal) nodes of UPGMA tree:")
    print(UPGMA_tree.get_nonterminals())


