# Dataframe-building related scripts and final taxonomy information spreadsheet

### Files included:
- Taxonomy information database.xlsx: All the species used in the datasets with Phylum, Order, family and genera specifie

- final_monophylo.txt: The same database, but compatable with MonoPhylo.py (Portik and Wiens, 2020)

### Scripts included:
- 01_Retrieve_classification_NCBI.py: returns the desired classification (phylum, order, family) based on inputted species name

- 02_Collect_sp_info.py: Collect species information from fasta files and store in a .csv file

- 03_Add_accessions.py: Where possible, retrieve the species accession numbers

- 04_In_concat_check.py: Return ```TRUE``` or ```FALSE``` depending on if the species exists within the 3-gene concatenation dataset (used for unconstrained and constrained analysis)





_Portik, D.M., and J.J. Wiens. (2020) Do alignment methods matter for phylogenomic (UCE) analyses? Systematic Biology, Advance Access_
