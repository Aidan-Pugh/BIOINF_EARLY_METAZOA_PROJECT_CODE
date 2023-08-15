# All files necessary to construct a constraint backbone tree, with the fianl constraint tree included

### Scripts:
- 01_Extract_orders.py: Searches the fasta file that will be used in the constrained analysis and creates a .txt file containing all the orders

- 02_Rename_newick_format.py: For each order in the order-level backbone tree, replace the order with all species belonging to that order and produce a new species-level constraint tree

### Files included:
- Backbone_newick_final.txt: Order-level constraint tree created in Mesquite and exported as a newick file
- Backbone_tree_species_level_revised.txt: Species-level constraint tree created from _02_Rename_newick_format.py_
- list_of_orders.txt: Output from _01_Extract_orders.py_
