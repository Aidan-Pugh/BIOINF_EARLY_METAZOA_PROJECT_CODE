# ps_test.fas

import os, argparse, re, csv
from Bio import Entrez, SeqIO, SeqRecord

Entrez.email = "mz22699@bristol.ac.uk"

######################################################
### FUNCTIONS FOR EXTRACTING INFORMATION FROM NCBI ###
######################################################

def get_tax_id(species):
    """ 
    Retrieve taxid from ncbi based on Species_name string
    
    This is useful for then calling the taxonomy information such 
    as order, using the function get_taxa_data below.
    """
    species = species.replace('_', "+").strip()
    
    search = Entrez.esearch(term = species, 
                            db = "taxonomy", 
                            retmode = "xml")
    record = Entrez.read(search)
    try:
        tax_id_number = record['IdList'][0]
        return tax_id_number
    except:
        print(f"No TaxID found for {species}")
        return None


def get_taxonomy_info(taxid):
    """
    Use the taxid name to fetch the record
    
    This is stored as an object containing class, order, 
    phylum, etc.
    """
    if taxid == None:
        return None
    else:
        try:
            result = Entrez.efetch(id = taxid, 
                                   db = "taxonomy", 
                                   retmode = "xml")
            tax_info_record = Entrez.read(result)
            return tax_info_record
        except:
            print(f"No taxonomy info found for {taxid}")
            return None


def extract_classification(tax_info_record, classification):
    """
    Extract the desired classification from the record
    
    Input is the output from "get_taxonomy_info"
    
    Takes the record producued by Entrez.efetch (function 
    "get_taxonomy_info") and returns the name of the desired 
    classification
    
    Inputs for the variable "classification" can include:
    superkingdom, clade, kingdom, phylum, class, subclass, order, 
    suborder, family, and genus
    
    """
    classification_name = "NA" # Automatically NA if nothing matches
    if tax_info_record == None:
        return classification_name
    else:
        for i in tax_info_record[0]["LineageEx"]:
            # Cycle through each rank until you find the classification specified
            if i["Rank"] == classification:
                classification_name = i["ScientificName"]
                print(f"The {classification} is {classification_name}")

        # If no classification is found, return NA for that classification
        if classification_name == "NA":
            print(f"No information found for classification {classification}")

        return classification_name


def get_the_order(species):
    """
    Combines the above functions into one package
    
    Accepts species, returns an order for that species
    """
    taxid = get_tax_id(species)
    tax_info = get_taxonomy_info(taxid)
    order=extract_classification(tax_info, "order")
    
    return order
    
    

def rename_headers(og_header, order, genus):
    new_header = og_header.replace(genus, order.upper()+"_"+genus)
    print(new_header)
    
    return new_header
    
    
    
    
# Function
def new_file(file_in, file_out_name, headers_list):
    file_name = open(file_in, "r")
    count = 0 
    newfasta = open(file_out_name, "w")
    
    for line in file_name:
        if line.startswith(">"):
            newname = headers_list[count]
            newfasta.write(newname+"\n")
            count += 1
        else:
            newfasta.write(line)
            

#############################################
### SCRIPT FOR UPDATING HEADER WITH ORDER ###
#############################################

print("Add orders to the headers")
#input_file = input("Please input the file name to update:")
genus = []
species = []
name = []
header = []
accession_nums = []


##############==> make lists of species, genus, headers,
with open("DONE_PORI.fas") as read:

    for line in read:
        if line.startswith(">"):
            ACSESSION = line.strip().lstrip(">")
            fields = ACSESSION.split("_")
            
            
            if len(fields) == 1:
                genus_only_temp = line.rstrip()
                genus_only = genus_only_temp.replace(">","")
                genus.append(genus_only)
                
                species_only = "sp"
                species.append(species_only)
                
                full_name_temp = "_".join([genus_only, species_only])
                
                name.append(full_name_temp)
                
                header.append(line.rstrip())
                
                
                
            if len(fields) == 2:
                genus_only_temp = line.split("_")[0]
                genus_only = genus_only_temp.replace(">","")
                genus.append(genus_only)
                
                species_only = "sp"
                species.append(species_only)
                
                full_name_temp = "_".join([genus_only, species_only])
                
                name.append(full_name_temp)
                
                header.append(line.rstrip())
                
                
            
            if len(fields) == 3:
            
                # Genus extraction
                genus_only_temp = line.split("_")[0]
                genus_only = genus_only_temp.replace(">","")
                genus.append(genus_only)

                # Species extraction
                species_only = line.split("_")[1]
                species.append(species_only)
                
                # Nums
                reference = line.split("_")[2]
                accession_nums.append(reference)
                
                # Full name
                full_name_temp = "_".join([genus_only, species_only])
                
                name.append(full_name_temp)

                # Header
                header.append(line.rstrip())
                
                
            if len(fields) >= 4:
                # Genus extraction
                genus_only_temp = line.split("_")[0]
                genus_only = genus_only_temp.replace(">","")
                genus.append(genus_only)
                # Species extraction
                species_only = line.split("_")[1]
                species.append(species_only)
                
                # Nums
                reference = line.split("_")[2]+"_"+line.split("_")[3]
                print(reference)
                accession_nums.append(reference)
                
                # Full name
                # Full name
                full_name_temp = "_".join([genus_only, species_only])
                name.append(full_name_temp)

                # Header
                header.append(line.rstrip())
            


############## Extract orders

order_list = []

for i in range(0,len(name)):
    order = get_the_order(name[i])
    order_list.append(order)

################# Rename headers with function:

new_headers_2 = []

for i in range(0, len(header)):
    renamed = rename_headers(header[i],order_list[i], genus[i])
    new_headers_2.append(renamed)


################## Write this all to a new file


new_file("DONE_PORI.fas", "DONE_PORI_ORDER.fas", new_headers_2)