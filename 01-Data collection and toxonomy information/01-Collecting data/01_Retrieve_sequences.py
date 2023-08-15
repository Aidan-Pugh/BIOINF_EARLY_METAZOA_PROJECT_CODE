#!/usr/bin/python
import os, argparse, re, csv
from Bio import Entrez, SeqIO, SeqRecord
# from Helpful_functions import *

Entrez.email = "mz22699@bristol.ac.uk"


######################################################
### FUNCTIONS FOR EXTRACTING INFORMATION FROM NCBI ###
######################################################

def get_tax_id(species):
    """ 
    Retrieve taxid from ncbi
    
    This is useful for then calling the taxonomy information such 
    as order, using the function get_taxa_data below.
    """
    species = species.replace('_', "+").strip()
    
    search = Entrez.esearch(term = species, 
                            db = "taxonomy", 
                            retmode = "xml")
    record = Entrez.read(search)
    
    return record['IdList'][0]


def get_taxonomy_info(taxid):
    """
    Use the taxid name to fetch the record
    
    This is stored as an object containing class, order, 
    phylum, etc.
    """
    result = Entrez.efetch(id = taxid, 
                           db = "taxonomy", 
                           retmode = "xml")
    tax_info_record = Entrez.read(result)
    
    return tax_info_record


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
    
    for i in tax_info_record[0]["LineageEx"]:
        # Cycle through each rank until you find the classification specified
        if i["Rank"] == classification:
            classification_name = i["ScientificName"]
            print(f"The {classification} is {classification_name}")
            
    # If no classification is found, return NA for that classification
    if classification_name == "NA":
        print(f"No information found for classification {classification}")
        
    return classification_name

##############################################
### FUNCTIONS FOR COLLECTING SEQUENCE DATA ###
##############################################

class TaxSeq:
    def __init__(self, identifier, recid, seq):
        self.identifier = identifier
        self.ids = [recid]
        self.seqs = [seq]
        self.taxid = None

    def __str__(self):
        out = ""
        name = self.identifier.split()

        out += f">{name[0]}_{name[1]}_{self.ids[0]}"
        out += f"\n{self.seqs[0]}"
        for i in range(len(self.ids[1:])):
            # ADDED THE ORDER
            out += f"\n>{name[0]}_{name[1]}_{self.ids[i]}\n"
            out += f"{self.seqs[i]}"

        return out

    def resolve_taxid(self):
        """ 
        Retrieve taxid from ncbi
        """
        search = Entrez.esearch(term = self.identifier, db = "taxonomy", retmode = "xml")
        record = Entrez.read(search)
        if len(record['IdList']) == 0:
            print("Failed to obtain taxid for " + str(self.identifier))
            return False
        self.taxid = record['IdList'][0]  
        print("Obtained tax id for " + str(self.identifier))
        return True  

# Get SeqRecords from NCBI Genbank
def get_entrez_sequences_iterator(searchTerm, retstart):
    """
    Get the sequences from NCBI
    """

    # Search for all the COI sequences in the family
    handle = Entrez.esearch(db="nucleotide", term=searchTerm, retstart=retstart*1500, retmax=1500)

    # Retrieve the list of matching sequence IDs
    record = Entrez.read(handle)
    id_list = record["IdList"]

    # Download the sequences in GB format
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")

    return records

def SeqRecords_iterator_to_TaxSeqs(recordsIterator, taxrecs):

    newseqs = 0
    totalrecs = 0

    for record in recordsIterator:
        totalrecs += 1

        seq = None
        # Extract COI feature(s)
        for feature in record.features:
            # Find features with CDS type and COI as gene name or product
            try:
                if feature.type == "CDS" and (feature.qualifiers["gene"][0].upper() in ["COI","COX1","CO1"] and feature.qualifiers["product"][0].lower() in ["cytochrome oxidase subunit 1","cytochrome c oxidase subunit i","cytochrome oxidase subunit i","cytochrome c oxidase subunit 1","cytochrome oxidase i"]):
                    # Extract the COI feature and create a new SeqRecord object
                    seq = str(feature.extract(record.seq))
                    
                    break

                # Extract organism name
                if "organism" in feature.qualifiers:
                    fullname = feature.qualifiers["organism"][0]
                    organism = fullname.split()[0] + " " + fullname.split()[1]

            except:
                continue
        
        # Extract taxon info if new entry

        if not seq:
            print("Record does not contain a COI feature. ID=" + record.id)
            continue

        if not organism:
            print("Record does not show organism source. ID=" + record.id)
        newseqs += 1

        if organism in taxrecs:
            taxrecs[organism].ids.append(str(record.id)) # Add the accession No
            taxrecs[organism].seqs.append(seq) # Add the sequence
        else:
            taxrecs[organism] = TaxSeq(organism, str(record.id), seq) # Create a new taxseq obj and supply info
    
    print(f"The total records looked at are {totalrecs}")

    return (taxrecs, newseqs)

# Collect for search term into taxrecs obj
def collect_and_parse_records(searchTerm):
    retstart = 0
    newseqs = 0
    taxrecs = {}
    # Collect and parse records
    while True:
        
        newseqs = 0
        recordsIterator = get_entrez_sequences_iterator(searchTerm, retstart)
        taxrecs, newseqs = SeqRecords_iterator_to_TaxSeqs(recordsIterator, taxrecs)

        print(f"Parsed {str(newseqs)} records from GenBank.")

        print(f"RUN: {retstart}. Parsed {newseqs} new records")

        retstart += 1

        if not newseqs: break
    
    return taxrecs

def new_file(taxrecs, filepath = None):
    if filepath == None:
        filepath = "out.fas"
    with open(filepath, 'w') as f:
        for taxid, rec in taxrecs.items():
            f.write(str(rec) + "\n")
    return True

def read_fasta_file_into_taxrecs(filepath):

    taxrecs = {}
    counter = 0

    with open(filepath, 'r') as f:
        # Define collection variables
        currentTaxID = ""
        currentSeqID = ""
        currentSeq = ""

        # Read lines in file
        for line in f:
            # If header
            if line.startswith(">"):
                # print("Header found")
                # If new header add last entry to list
                if currentTaxID:
                    counter += 1 
                    if currentTaxID in taxrecs:
                        taxrecs[currentTaxID].ids.append(currentSeqID)
                        taxrecs[currentTaxID].seqs.append(currentSeq)
                    else:
                        # print("New Taxa")
                        taxrecs[currentTaxID] = TaxSeq(currentTaxID, currentSeqID, currentSeq)
                
                currentSeq = ""                             
                #Get header
                header = line.strip().lstrip(">")
                fields = header.split("_")
                # print(header)
                # print(fields)
                
                #Collect details
                if len(fields) == 3:
                    currentTaxID = fields[0] + ' ' + fields[1]
                    currentSeqID = fields[2]
                elif len(fields) == 4:
                    currentTaxID = fields[0] + ' ' + fields[1]
                    currentSeqID = fields[2] + "_" + fields[3]
                elif len(fields) == 2:
                    currentTaxID = fields[0] + ' ' + fields[1]
                    currentSeqID = "N"
                else:
                    currentTaxID=""

            else:
                #Read sequence
                currentSeq += line.strip()
        
        # add the last entry
        if currentTaxID:
            counter += 1  
            if currentTaxID in taxrecs:
                taxrecs[currentTaxID].ids.append(currentSeqID)
                taxrecs[currentTaxID].seqs.append(currentSeq)
            else:
                # print("New Taxa")
                taxrecs[currentTaxID] = TaxSeq(currentTaxID, currentSeqID, currentSeq)

    print(f"Successfully parsed {str(counter)} records from {filepath}.")

    # print(taxrecs)

    return taxrecs

def prune_by_largest_coverage(taxrecs, min, max, optima):

    for organism, rec in taxrecs.items():
        
        # Pruning by hard limits where min=XXXX, and max = XXXX
        pruned = []
        for i in range(len(rec.ids)):
            if len(rec.seqs[i]) < max and len(rec.seqs[i]) > min:
                pruned.append([rec.ids[i], rec.seqs[i]])

        if not pruned:
            print(f"No sequences within the hard limits could be obtained for {organism}. \nPlease inspect manually.")
            for i in range(len(rec.ids)):
                pruned.append([rec.ids[i], rec.seqs[i]])

        # Find sequence closest to optima
        closest = None
        closestDist = float('inf')
        for i in range(len(pruned)):
            # Looks at modulus of seq length minus optimal length
            distance = abs(len(pruned[i][1]) - optima)
            if distance < closestDist:
                closest = pruned[i]
                closestDist = distance

        # We are left with a tuple (id, seq) which is closest to the optimal length
        # Set ids and seqs to nothing; replacing with our closest id seq pair

        rec.ids = [closest[0]]
        rec.seqs = [closest[1]]

    return taxrecs
        

###############################################
#### SCRIPT FOR THE RETRIEVAL AND RENAMING ####
###############################################


Entrez.email = "mz22699@bristol.ac.uk"
    
# "Parazoanthus swiftii"[Organism] AND coi[All Fields] AND "complete genome" 
# "Acropora abrolhosensis"[Organism] AND (COI OR CO1 OR Cytochrome C oxidase)
# Heteropolypus
# "Alcyoniidae"[Organism] AND COI[All Fields] 

# Lines to collect records:
#print("Running function 'collect_and_parse_records'")
# taxrecs = collect_and_parse_records('"Cnidaria"[Organism] AND (COI[gene] OR cox1[gene] OR co1[gene] OR coxI[gene])')

# print("Collection completed")

#Display stats
#totalnoseqs = 0
#for organism, rec in taxrecs.items():
    # print(f"{organism} : {len(rec.ids)} sequences") # Uncomment for per taxa stats DO NOT USE WITH 1000s OF TAXA
    #totalnoseqs += len(rec.ids)

#print(f"Found {str(totalnoseqs)} over {str(len(taxrecs.keys()))}")


# Lines to prune records and save to file


print("Reading file")
taxrecs = read_fasta_file_into_taxrecs("CNID_COI_RAW.fas")
print("Pruning")
taxrecs = prune_by_largest_coverage(taxrecs, 100, 1400, 750)
print("Writing the results to file")
new_file(taxrecs, "CNID_COI_PRUNED.fas")

#Display stats
totalnoseqs = 0
for organism, rec in taxrecs.items():
    # print(f"{organism} : {len(rec.ids)} sequences") # Uncomment for per taxa stats DO NOT USE WITH 1000s OF TAXA
    totalnoseqs += len(rec.ids)

print(f"Found {str(totalnoseqs)} over {str(len(taxrecs.keys()))}")

# Line for printing to a file
#print("Writing the results to file")
#new_file(taxrecs, "CNID_COI_RAW.fas")

