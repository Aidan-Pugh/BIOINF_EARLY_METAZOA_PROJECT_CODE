# Collecting data from NCBI directly using biopython and Entrez

For collecting COI data:
```if feature.type == "CDS" and (feature.qualifiers["gene"][0].upper() in ["COI","COX1","CO1"] and feature.qualifiers["product"][0].lower() in ["cytochrome oxidase subunit 1","cytochrome c oxidase subunit i","cytochrome oxidase subunit i","cytochrome c oxidase subunit 1","cytochrome oxidase i"]):```

For collecting 18S data:```if feature.type == "rRNA" and (feature.qualifiers["product"][0].lower() in ["18s ribosomal rna"]):```

For collecting CYTB data:```if feature.type == "CDS" and (feature.qualifiers["gene"][0].upper() in ["CYTB"] and feature.qualifiers["product"][0].lower() in ["cytochrome b"]):```
