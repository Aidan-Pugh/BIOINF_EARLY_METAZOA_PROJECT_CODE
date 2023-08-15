import pandas as pd
import csv

with open("DATAFRAME.csv","r") as taxa_file:
    reader = csv.reader(taxa_file)

    species_list = []

    for row in reader:
        # row with species information == 4
        species_list.append(row[4])

# remove the first "species" as this is the column header
species_list.pop(0)


concat_species = []
with open("concat_data.fas") as concat:
    for line in concat:
        if line.startswith(">"):
            fields = line.split("_")
            species_name = "_".join([fields[1], fields[2]]).rstrip("\n")
            concat_species.append(species_name)

constraint_tree_list = []

# for each species in the dataframe:
for s in species_list:
    # if species in the concat data:
    if s in concat_species:
        tf = "TRUE"
    else:
        tf = "FALSE"
    constraint_tree_list.append(tf)

# double check any issues where a concat species is not in the database
# this shouldn't print anything
for s in concat_species:
    if s not in species_list:
        print(s, "not in species list")


df = pd.read_csv("DATAFRAME.csv")
# Add TRUe and FALSE to the dataset
df['in_constraint_tree'] = constraint_tree_list

# Create new .csv file with concat information in
df.to_csv("OUTPUT.csv")