#!/usr/bin/env python
import re 

file1 = "Homo_sapiens.GRCh38.pep.all.fa" #open Homo Sapiens protein fasta file downloaded from Ensembl.org
#file1 = "Danio_rerio.GRCz11.pep.all.fa" #open Danio Rerio protein fasta file downloaded from Ensembl.org

gene_id_length = {} #create dictionary to hold gene name as key and longest associated protein sequence as value 
gene_headers = {} #create a dictionary to hold new fasta headers for output file 

with open(file1, "r") as fh:   #open file 
    
    lines = fh.readlines()  #store lines of file in memory 
    
    for i, line in enumerate(lines):
        
        if line.startswith('>'):   #if the current line is a header line
            gene = line.split(" ")[3].split(":")[1].split(".")[0]    #split the contents of the line, retain the gene ID in a variable 
            protein = line.split(" ")[0].split(".")[0]               #retain the protein ID in a variable 
            gene_match = re.search("gene_symbol:[^ ]+", line)        #search for gene name 
            
            if gene_match != None: #if gene name found 
                gene_name = gene_match.group(0).split(':')[1]    #strip gene match of unecessary charachters and assign to gene name variable 
            else:
                gene_name = "No Name"     #if gene name not found, assign "No name" to gene name variable  
            
            header = protein + " " + gene + " " + gene_name   #create header with protein ID gene ID and gene name 
            
            i += 1 

            sequence = ""   #variable to store sequences while searching for the longest per gene ID 
            while (i < len(lines)) and (not lines[i].startswith('>')):  #while not at the end of the file and not on a header line
                sequence = sequence + lines[i].strip('\n')    #append line to sequence variable 
                i += 1
            
            if gene not in gene_id_length and gene not in gene_headers:    #if geneID doesnt exist in any of the dictionaries 
                gene_id_length[gene] = sequence     #add it as a new entry to both 
                gene_headers[gene] = header
            elif gene in gene_id_length and gene in gene_headers:   #if gene ID exists in the dictionaries 
                if len(gene_id_length[gene]) < len(sequence): #and if the length of the corresponding sequence value of the geneID key is less than the length of the current sequence 
                    gene_id_length[gene] = sequence #replace the current sequence with the old one
                    gene_headers[gene] = header  #update the header in the header dictionary to reflect this change 

    f = open("homo_sapien_longest_prot.fa", "w")   #open output file 

    for key in gene_headers:    #write new headers and corresponding longest protein sequences per gene to the file 
        f.write(gene_headers[key])
        f.write("\n")
        if key in gene_id_length:
            f.write(gene_id_length[key])
            f.write( "\n")

    f.close()
