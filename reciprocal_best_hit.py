import os 
from os import path
import gzip

bestevalue = {}
rbh = {}
prev_geneID = ''
prev_evalue = 0

output = open('Human_Zebrafish_RBH.tsv', 'w')
hg = open('human_genes.txt', 'r')
zg = open('zfish_biomart.txt', 'r')

with open('zfish_blastp_human_sorted.txt') as f2h:

    line = f2h.readline().strip('\n') 
    
    while 1:
        line = f2h.readline().strip('\n')
        if not line: 
            break
        geneIDF = line.split()[0] 
        gene_IDH = line.split()[1]
        e_value = line.split()[10]
        pair = gene_IDH + geneIDF
        
        if geneIDF != prev_geneID:
            bestevalue[pair] = line
        
        if geneIDF == prev_geneID:
            if e_value == prev_evalue:
                bestevalue.pop(pair, None)
            else:
                continue

        prev_geneID = geneIDF
        prev_evalue = e_value

with open('human_blastp_zfish_sorted.txt') as h2f:

    line = h2f.readline().strip('\n') 
    
    while 1:
        line = h2f.readline().strip('\n')
        if not line: 
            break
        geneIDH = line.split()[0]
        geneIDF = line.split()[1]
        e_value = line.split()[10]
        pair = geneIDH + geneIDF

        if geneIDH != prev_geneID:
            if pair in bestevalue:
                if e_value == bestevalue[pair].split('\t')[10]:
                    rbh[bestevalue[pair]] = line 
                else:
                    bestevalue.pop(pair, None) 

        if geneIDH == prev_geneID:
            if pair in bestevalue:
                if e_value == prev_evalue:
                    rbh.pop(bestevalue[pair], None)
                else:
                    continue        

        prev_geneID = geneIDH
        prev_evalue = e_value

    for key in rbh:
        print("Fish: " + key)
        print("Human:" + rbh[key])


output.write("Human Gene ID\tHuman Protein ID\tHuman Gene Name\tZebrafish Gene ID\tZebrafish Protein ID\tZebrafish Gene Name\n")

for key in rbh:
    fish_prot = str(key).split()[0]
    human_prot = str(rbh[key]).split()[0]

    hg.seek(0)
    zg.seek(0)

    for line in hg:
        if len(line.split()) == 3 and line.split()[1] == human_prot: 
            print(len(line.split()))
            print(human_prot)
            for line2 in zg:
                if len(line2.split()) == 3 and line2.split()[1].strip('\t') == fish_prot:
                    output.write(line.strip('\n') + "\t" + line2)

output.close()
