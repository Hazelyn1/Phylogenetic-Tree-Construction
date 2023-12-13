# Final_project
This program constructs a UPGMA tree and NJ tree using the BioPython package. This program takes the results of a DNA multiple sequence alignment (via online tool MUSCLE) in .clw format as input, which is provided in an example file. The distance matrix is calculated first and then the "four_point_condition" function is called which checks if the distance matrix is additive or not by checking the following formula:
max{𝑑(𝑎, 𝑐) + 𝑑(𝑏, 𝑑), 𝑑(𝑎, 𝑑) + 𝑑(𝑏, 𝑐)} ≥ 𝑑(𝑎, 𝑏) + 𝑑(𝑑, 𝑐)
The results of the "four_point_condition" function do not impact tree construction since this program implements heuristic methods. 
The trees get constructed by creating an object of the type DistanceTreeConstructor and drawn to a new window using the ".draw" function. The program prompts the user which tree they would like to see and also provides the terminal and non-terminal branch lengths of the tree(s) to the user.
