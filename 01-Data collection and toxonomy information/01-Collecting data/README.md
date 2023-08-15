# Collecting data from NCBI directly using biopython and Entrez
### 01_Retrieve_sequences.py
- Requires updating the species/ phylum of interest and the gene query.
- For collecting _COI_ data:
```if feature.type == "CDS" and (feature.qualifiers["gene"][0].upper() in ["COI","COX1","CO1"] and feature.qualifiers["product"][0].lower() in ["cytochrome oxidase subunit 1","cytochrome c oxidase subunit i","cytochrome oxidase subunit i","cytochrome c oxidase subunit 1","cytochrome oxidase i"]):```
- For collecting _18S_ data:```if feature.type == "rRNA" and (feature.qualifiers["product"][0].lower() in ["18s ribosomal rna"]):```
- For collecting _CYTB_ data:```if feature.type == "CDS" and (feature.qualifiers["gene"][0].upper() in ["CYTB"] and feature.qualifiers["product"][0].lower() in ["cytochrome b"]):```
- Saves everything to a .txt file in fasta format. Uncomment section of the script to prune to one sequence per species

### 02_Append_order.py
- Search the fasta file for species names. Use species names to retrieve order classficaiton information from NCBI taxonomy
- Append order to the start of the headers to make future analysis easier
- Saved to a new file
