import re
import os

######################################################
### 1. FUNCTIONS - DO NOT EDIT ANY OF THIS SECTION ###
######################################################

def get_classification_genus_species(line_in, classification_dict):
    """
    Take in a line, extract the classification and the species, and add this info to a dictionary.
    
    E.g: get_classification_genus_species(">SCINCIDAE_Ablepharus_bivittatus", family_dictionary)
         print(family_dictionary)
         {"SCINCIDAE" : "Ablepharus_bivittatus"}
         
    """
    # Pattern of the fasta file
    pattern = r"^>([A-Z]+)_(\w+)$" # Assumes the following: >FAMILY_Species_name
    match = re.match(pattern, line_in)
    
    # extract the classification and the species name
    if match:
        classification = match.group(1)
        genus_species = match.group(2)
        classification_genus_species = "_".join([classification, genus_species])
        
        # Add the classification as a key and the species name as values
        if classification not in classification_dict:
            classification_dict[classification] = []
        if not any(classification_genus_species in values for values in classification_dict.values()):
            classification_dict[classification].append(classification_genus_species)
            
def fill_dictionary(fasta_file, classification_dict):
    """
    Look through the fasta_file and extract all the species names into collection of families
    
    E.g: {"FamilyA" : ["Species_a", "Species_b", "Species_c"],
          "FamilyB" : ["Species_d", "Species_e", "Species_f"]...}
    """
    with open(fasta_file) as read:
        for line in read:
            get_classification_genus_species(line, classification_dict)
    return print(f">>>> Info collected in dictionary for file {fasta_file}<<<<")  

def update_and_write(tree_file, classification_dict_full):
    """
    tree_file is the name of newickfile you are editing
    classification_dict_full is the dictionary where family and species info was stored
    
    This function writes everything to a new file
    """
    with open(tree_file) as file:
        text = file.read()
        for key in classification_dict_full:
            if key in text:
                values = classification_dict_full.get(key)
                new_values = ",".join(values)
                updated = text.replace(key, new_values)
                text = updated
    return updated



##################################################
### 2. FILE NAMES - SOME EDITS ARE NEEDED HERE ###
##################################################

# This is the name of your backbone produced by mesquite that you are editing
tree_file = "Newick_backbone.txt"

#This is the name of your final and updated newick format - change if you want but this is fine as is
new_file = "0_Backbone_tree_species_level.txt"

# This is the concatenated data with the families appended on the front
fasta_file = "concat_data.fas"


###################################################
### 3. SCRIPT - DO NOT EDIT ANY OF THIS SECTION ###
###################################################

family_dict = {}

fill_dictionary(fasta_file, family_dict)
print(len(family_dict))
print(sum(len(values) for values in family_dict.values()))

with open(new_file, "w") as file:
    file.write(update_and_write(tree_file, family_dict))
    
print(">>>> Written to a new file <<<<")
print(f">>>> Load the new_file into FigTree to make sure it's correct <<<<")