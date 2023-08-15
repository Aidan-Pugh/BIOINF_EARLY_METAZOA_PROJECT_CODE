#!/bin/bash

# Run IqTree for the 3-gene concatenated data, specifiying a file to be used as a constraint tree

iqtree2 -s FcC_supermatrix_18S_CYTB_COI_Constrained.phy -g Backbone_tree_species_level_revised.txt -m MFP -B 1000 -nt AUTO -mem 20G

exit