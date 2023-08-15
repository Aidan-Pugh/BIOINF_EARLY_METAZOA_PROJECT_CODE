## THIS WORKS PERFECTLY!
# NOW, TURN IT INTO A FUNCTION SO I CAN RUN IT FOR ALL THREE GENES

import pandas as pd

def read_gene_tree_files(file_in, list_species, list_accession):
    """
    Function to read through fasta file and extract the accession number
    for the species
    """
    with open(file_in) as read:
        for line in read:
            if line.startswith(">"):
                fields = line.split("_")
                if len(fields) >= 3:
                    if fields[2] == fields[-1]:
                        genus_species = "_".join([fields[1], fields[2].rstrip("\n")])
                        accession = "CHECK"
                        list_species.append(genus_species)
                        list_accession.append(accession)
                    else:
                        genus_species = "_".join([fields[1], fields[2]]).rstrip("\n")
                        accession = fields[-1].rstrip("\n")
                        list_species.append(genus_species)
                        list_accession.append(accession)
                else:
                    genus_species = "_".join([fields[1], "SPECIES", "CHECK"]).rstrip("\n")
                    list_species.append(genus_species)
                    list_accession.append(accession)
            

    return print(f"{file_in} read complete")


# Read in all species from taxa dataframe
df = pd.read_csv("TEST_COI_18S.csv", header = 0)
print(df.iloc[0:10])
species_list_all = df["Species"].tolist()

# store species list for cytb & accessions for 18s
gene_tree_18s_list = []
accessions_list_18s = []
read_gene_tree_files("CYTB_gene_tree.fas", gene_tree_18s_list, accessions_list_18s)

# Store the final list of accession for the all species dataset
Final_accession_col = []
success = []

for s in species_list_all:
    if s in gene_tree_18s_list:
        success.append(s)
        index = gene_tree_18s_list.index(s)
        accession_num = accessions_list_18s[index]
        Final_accession_col.append(accession_num)
    else:
        Final_accession_col.append("FALSE")


df["CYTB accession"] = Final_accession_col

df.to_csv("TEST_COI_18S_CYTB.csv")