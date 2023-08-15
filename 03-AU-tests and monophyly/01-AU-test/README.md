# Directory containing the scripts used to prune the trees and run the AU-tests, as well as the outputted results

### Scripts:
- 01_Run_AUtest.sh: Script that calls the other two scripts. Loops through all directories starting ```C_C_``` to complete the AU tests all in one sitting
  
- 02_Prune_trees.py: Python script to prune tree in the same folder to shared taxa, a prerequisite to AU-tests being run

- 03_AUtest.sh: Script for the runnning of the AU-test

### Other files:
- File_structure.png: Image of the required directory layout
  
- AU test results for table.xlsx: Raw data produced from some of the AU tests

- Sample_output.zip: Example of the output after running the above script on one directory of trees


### File structure when running script: IMPORTANT
- The only script to be run is 01_Run_AUtest.sh, from a directory containing directories of interest (all must start with C_C_)
- Within each tree directory, the following must be present:
  - 3-gene unconstrained file (.contree) named ```02_concat_tree.phy.contree```
  - Tree(s) of interest (.contree)
  - Fasta file used to construct the unconstrained tree
