# Files used to construct the backbone tree for the constrained analysis

1. extract_order.py used to extract a list of the orders from a stated fasta file that you will be constraining
2. Rename_newick.py used to fill in the orders with species of that order, based on the fasta file being constrained

Other files included:
- list_of_orders.txt                     Orders extracted using the extract_order.py script
- backbone_mesquite_format               Backbone tree in the format exported from Mesquite
- Backbone_newick_final.txt              The order-level constraint tree created in Mesquite and exported in Newick format
- Backbone_tree_species_level_revised    The species-level constraint tree filled in using Rename_newick.py
