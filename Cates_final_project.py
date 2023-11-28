#Hazelyn Cates
#Started 11/6/23, finished 11/16/23
#EN.605.651.81.FA23
#This program constructs phylogenetic trees using Biopython

import Bio
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

#function that checks if the calculated distance matrix is additive
def four_point_condition(dist_matrix):
    from itertools import permutations
    import math

    flag = 0 #flag to keep track of if matrix is additive. 0 = not additive, 1 = additive
    #First, calculate how many combinations of 4 there are:
    print("Enter number of sequences:")
    n = int(input())
    k = 4

    num_combos = math.comb(n, k)
    #print(num_combos)

    names = dist_matrix.names #accession numbers

    #Now find all the accession number combos of 4:
    perms = permutations(names, 4)
    combos = [] #array to hold the permutations

    for i in list(perms): #add all permutations to a list to create a list of lists
        combos.append(i)

    #Inequaility to check if matrix is additive: max{𝑑(𝑎, 𝑐) + 𝑑(𝑏, 𝑑), 𝑑(𝑎, 𝑑) + 𝑑(𝑏, 𝑐)} ≥ 𝑑(𝑎, 𝑏) + 𝑑(𝑑, 𝑐)
    for i in range(0, num_combos): #number of possible permutations of the 10 sequences (210 total)
        a = names.index(combos[i][0]) #i = index of permutation sequence with the "combos" list (it's a list of lists)
        b = names.index(combos[i][1])
        c = names.index(combos[i][2])
        d = names.index(combos[i][3])

        val_1 = max((round(dist_matrix[a][c], 1) + round(dist_matrix[b][d], 1), round(dist_matrix[a][d],1) + round(dist_matrix[b][c], 1)))
        #print(val_1)
        val_2 = round(dist_matrix[a][b], 1) + round(dist_matrix[d][c], 1)
        #print(val_2)
        if val_1 >= val_2: #So far, the matrix is additive
            flag = 1
            continue

        else: #matrix is not additive
            flag = 0
            break #end for loop

    #check state of flag
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
#To do this, first calculate the distances between the sequences:
calc = DistanceCalculator('identity') #for DNA distance calculations

#Now generate the distance matrix of the aligned sequences using this distance calculator
dist_matrix = calc.get_distance(aligned_seqs)

#Determine if distance matrix is additive
four_point_condition(dist_matrix)

#From here, the NJ and UPGMA algorithms can be used to generate trees and compare the results
construct = DistanceTreeConstructor()

NJ_tree = construct.nj(dist_matrix)

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
        print("\nBranch lengths of terminal nodes (leaves) of UPGMA tree:")
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
        print("\nBranch lengths of terminal nodes (leaves) of NJ tree:")
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
    print("\nBranch lengths of terminal nodes (leaves) of UPGMA tree:")
    print(UPGMA_tree.get_terminals())
    print("Branch lengths of non-terminal (internal) nodes of UPGMA tree:")
    print(UPGMA_tree.get_nonterminals())


