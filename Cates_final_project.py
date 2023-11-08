#Hazelyn Cates
#Started 11/6/23
#EN.605.651.81.FA23
#This program constructs a phylogenetic tree using Biopython package

#DELETE THE FOLLOWING COMMENT, ADD TO REFERENCES IN REPORT!!!!
#Largely based on this tutorial: https://medium.com/@poudelmohit59/beginners-guide-to-phylogenetic-tree-construction-using-biopython-5accbd8345a2

import math
import Bio
from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def four_point_condition(dist_matrix): #checks if the calculated distance matrix is additive
    from itertools import permutations
    #First, calculate how many combinations of 4 there are:
    print("Enter number of sequences:")
    n = int(input())
    k = 4

    num_combos = math.comb(n, k) #https://www.w3schools.com/python/ref_math_comb.asp
    #print(num_combos)

    #Now find all the combos of 4: https://www.geeksforgeeks.org/permutation-and-combination-in-python/
    names = dist_matrix.names
    perms = permutations(names, 4)
    combos = [] #array to hold the 210 permutations

    for i in list(perms): #add all permutations to a list
        print(i)
        combos.append(i)
    #print(combos)

    for i in range(0, len(combos)): #number of possible permutations of the 10 sequences
        for j in range(0, 4): #number of sequences in each permutation
            print(combos[7][2]) #i = index of permutation sequence, j = index of specific accession # in that permutation
            #For for example you have the first permutation: ('NM_001081819.2', 'NM_000594.4', 'NM_001003244.4', 'NM_214022.1')
            #Which is at index i = 0
            #Say you want to access the third accession # in that permutation sequence. It would be j = 2
            #So it's like a list of lists

    #Now I need to apply the following formula
    #max{𝑑(𝑎,𝑐)+𝑑(𝑏,𝑑),𝑑(𝑎,𝑑)+𝑑(𝑏,𝑐)}≥𝑑(𝑎,𝑏)+𝑑(𝑑,𝑐)
    #Where a = 1, b = 2, c = 3, d = 4
    #So when accessing these combination, need to access a positin WITHIN the list (is that possible?? idk)
    #So like in the distance matrix, I need to access the correct values given the accession numbers in the current combo of 4
    print(dist_matrix)
    #print(dist_matrix[3, 2])

    #for example if you have the combo: ('NM_001081819.2', 'NM_000594.4', 'NM_001003244.4', 'NM_214022.1')
    #Then you have a, b, c, and d
    #So the distance between a (row 9) and c (column 7) would be at index [9, 7]


    #So I image I would need a loop go through through all combos
    #But I somehow need to assign a, b, c, and d in each combo and each iteration
    max_vals = [] #empty array to store the "max" part of the above formula
    comparison_vals = [] #empty array to store the "greater than or equal to" part of the formula

    #Both arrays should have 210 values
    #Which will allow for the comparison
    #And if just one value of the 210 does NOT satify the forumla, the matrix is NOT additive and the loop can stop

    """for i in range(0, num_combos): #iterates 210 times
        print(list(combos[i]))"""





#read in FASTA sequences from file ad put them in an array that can be indexed:

#10 FASTA sequences were aligned via the online tool MUSCLE (https://www.ebi.ac.uk/Tools/msa/muscle/)
#The alignment file (in .clw format) is read in using AlignIO package:
with open("msa_muscle_output.clw", "r") as msa:
    aligned_seqs = AlignIO.read(msa, "clustal")

#print(aligned_seqs)

#FIX THIS!!! Need to somehow change the names of each sequence b/c it's just the accession numbers rn and I have no idea which is which
#But if I can't, I'll make sure I specify what's what in my report when I do my results section

#From the alignment data, have to create a distance matrix
#To do this, first calculate the distances between the sequences: (NOTE: I'm at part 3 in the tutorial linked above)
calc = DistanceCalculator('identity') #for DNA distance calculations

#Now generate the distance matrix of the aligned sequences using this distance calculator
dist_matrix = calc.get_distance(aligned_seqs)
four_point_condition(dist_matrix)
#print(dist_matrix)

#FIX THIS!!!! NEED TO ADD SOMETHING THAT CHECKS IF THE MATRIX SATISFIES THE FOUR-POINT-CONDITION!! need to make sure it's additive
#If it is, continue on.
#If not, end the program and ask user for new input



#From here, the NJ and UPGMA algorithms can be used to generate trees and compare the results
#First, the NJ method: https://homolog.us/Biopython/Bio.Phylo.TreeConstruction.DistanceTreeConstructor.html#content
construct = DistanceTreeConstructor()
NJ_tree = construct.nj(dist_matrix)
#Phylo.draw(NJ_tree)

#Repeat for UPGMA method:
UPGMA_tree = construct.upgma(dist_matrix)

#print(UPGMA_tree)
#Phylo.draw(UPGMA_tree)




