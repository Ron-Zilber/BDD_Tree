# BDD_Tree
The provided Erlang code builds a reduced Binary Decision Diagram (BDD) tree (https://en.wikipedia.org/wiki/Binary_decision_diagram) from a given boolean function. The returened tree will be the most efficient tree, by an ordering criterion input from the user.
At first, all possible permutations of trees are built in map or record data structure (by the user's choice).
Afterwards, the most efficient tree is found and returned to the user. 
The most efficient tree will be chosen by one of three possible user inputs: num_of_leaves, num_of_nodes or tree_height.

This project was writen as a middle assignment in an acadamic course -  'Functional Programming in Concurrent and Distributed Systems' for Electrical and Computers Engineering students at Ben-Gurion university. 
The requirements are detailed in the .pdf file

# Usage:
The code supports four functions:
1. exp_to_bdd(BoolFunc, Ordering, DataStructure) -> BddTree:
  this functions accepts a tuple BoolFunc by the structure that defined in the pdf file,
  desired Ordering criterion: tree_height, num_of_nodes or num_of_leafs (atoms) and DataStructure - map or record (atoms), and returns the most efficient tree.

2. solve_bdd(Bdd_tree, List) -> Result:
  This function accepts a BDD tree that has been return from exp_to_bdd function, and a list fo values that correspond to the variables. The list contains tuples in the form of [{x1, val1}, {x2, val2},..., {xn_valn}] where the values are 0, 1, false or true.
The result will be false or true.

3. listOfLeaves(Bdd_tree) -> Res:
  This function accepts a BDD tree and returns a list of lists, such that each list is a path from the root to a leaf. (Each leaf produce such list).

4. reverseIteration(LeafPtr) -> Res:
  This function accepts a list ("pointer") that has been returened from listOfLeaves function, and returns the list of nodes on the shortest path to the root.
  
  
  to use this code, compile the BDD_Tree.erl file on your machine, create a boolean function by the structure that has been defined in the pdf file and call function (1) from you linux shell. After this step, use the returened tree as an input to function (2), (3), and use the returned value from function (3) as an input to function (4).
