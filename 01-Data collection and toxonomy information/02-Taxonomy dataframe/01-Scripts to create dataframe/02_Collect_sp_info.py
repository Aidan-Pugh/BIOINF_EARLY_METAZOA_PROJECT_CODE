import pandas as pd

def read_gene_tree_files(file_in, header_list):
    """
    Read a file and add the order, genus and species to a list for each header
    """
    with open(file_in) as read:
        for line in read:
            if line.startswith(">"):
                line.rstrip(">")
                fields = line.split("_")
                if len(fields) >= 3:
                    order_genus_species = "_".join([fields[0],fields[1],fields[2]]).rstrip("\n")
                    if order_genus_species not in header_list:
                        header_list.append(order_genus_species)
                    else:
                        print(f"{order_genus_species} already in the list")
                else:
                    order_genus_species = "_".join([fields[0], fields[1], "SPECIES"]).rstrip("\n")
    return print(f"{file_in} read complete")

##########
# SCRIPT #
##########

header_list = []
# All files to read:
all_species_files = ["18s_gene_tree.fas", "COI_gene_tree.fas","CYTB_gene_tree.fas"]

for f in all_species_files:
    read_gene_tree_files(f, header_list)

Order = [] # MALACALCYONACEA
Genus = [] # Incrustatus
Species = [] # Incrustatus_comauensis
Accession = [] #KW13242

for h in header_list:
    fields = h.split("_")
    # Get the order
    order_temp = fields[0].lstrip(">")
    Order.append(order_temp)
    print(order_temp)

    # Get the genus
    genus_temp = fields[1]
    Genus.append(genus_temp)
    print(genus_temp)

    # get the species
    species_temp = fields[1] + "_" + fields[2]
    Species.append(species_temp)
    print(species_temp)

# Create a database to add this information to
taxon_dataframe = pd.DataFrame()

# Add the classifications
taxon_dataframe["Order"] = Order
taxon_dataframe["Genus"] = Genus
taxon_dataframe["Species"] = Species

taxon_dataframe.to_csv("TAXA_FILE_OUTPUT.csv")