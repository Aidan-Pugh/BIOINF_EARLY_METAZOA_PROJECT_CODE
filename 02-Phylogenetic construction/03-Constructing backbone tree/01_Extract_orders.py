### THIS SCRIPT EXTRACTS THE order FROM THE FASTA FILE(S) 
### EACH order WRITTEN TO A FILE ONCE IN BRACKETS
### Rename_newick_format.py CAN THEN BE USED TO CRAETE A BACKBONE

import re
import os

def extract_order(line_in, order_list):
    """
    Retrieve the order and add to a list
    
    No repeated orders in the list
    """
    pattern = r"^>([A-Z]+)_(\w+)$"
    match = re.match(pattern, line_in)
    
    if match:
        order = match.group(1)
        order_brackets=f"({order})"
        if order_brackets not in order_list:
            order_list.append(order_brackets)


def read_file(file_name, order_list):
    """
    Open the file to be read, runs extract_order
    """
    with open(file_name) as read:
        for line in read:
            extract_order(line, order_list)
    
def write_to_file(file_name, list_to_write):
    """
    Write the orders list to .txt file
    """
    with open(file_name, "a") as write_file:
        write_file.write("(")
        write_string = ",".join(list_to_write)
        write_file.write(write_string)
        write_file.write(")")
    
def run(class_list, file_read, file_write):
    """
    Read and write orders to new file
    """
    read_file(file_read, class_list)
    write_to_file(file_write, class_list)
    return print(">>>> RUN")

file_list = ["please_work.fas"]

order_list = []
for f in file_list:
    print(f">>>> Reading {f}")
    read_file(f, order_list)
write_to_file("All_families.txt", order_list)
print(f"Number of orders written:{len(order_list)}")