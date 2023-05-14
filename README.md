# BDD_Tree
The provided Erlang code builds a reduced Binary Decision Diagram (BDD) tree (https://en.wikipedia.org/wiki/Binary_decision_diagram) from a given boolean function. It recursively constructs the BDD by splitting the tree based on the variables in the function. 
At first, all possible permutations of trees are built in map or record data structure (by the user's choice), and afterwards, the most efficient tree are found and returned to the user. The most efficient tree will be chosen by one of three possible user inputs: num_of_leaves, num_of_nodes or tree_height.
This project was writen as a middle assigment in an acadamic course 'Functional Programming in Concurrent and Distributed Systems' for Electrical and Computers Engineering students at Ben-Gurion university. 
The requirements are detailed in the .pdf file
