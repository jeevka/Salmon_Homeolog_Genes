import sys
import re

F1 = open("Salmon_Homeologs_051214.gff3","r")


CDS  = {}
mRNA = {}
UTR5 = {}
UTR3 = {}

for i in F1:
    temp = i.split()
    if temp[2] == "CDS":
        ID = temp[8].replace("Parent=","")
        CDS[ID] = ID
       
    if temp[2] == "five_prime_UTR":
        ID = temp[8].replace("Parent=","")
        UTR5[ID] = ID
        
    if temp[2] == "three_prime_UTR":
        ID = temp[8].replace("Parent=","")
        UTR3[ID] = ID
        
    if temp[2] == "Homeolog_Transcript":
        temp2 = temp[8].split(";")
        ID = temp2[0].replace("ID=","")
        mRNA[ID] = ID

F1.seek(0)

# IDS to delete
IDS = []
for i in CDS:
    if not mRNA.has_key(i):
       IDS.append(i)

for i in UTR5:
    if not mRNA.has_key(i):
       IDS.append(i)

for i in UTR3:
    if not mRNA.has_key(i):
       IDS.append(i)


for i in F1:
    temp = i.split()
    
    if temp[2] == "CDS":
        ID = temp[8].replace("Parent=","")
        CDS[ID] = ID
    
    if temp[2] == "Homeolog_Transcript":
        temp2 = temp[8].split(";")
        ID = temp2[0].replace("ID=","")
        
    if temp[2] == "three_prime_UTR":
        ID = temp[8].replace("Parent=","")
        
    if temp[2] == "five_prime_UTR":
        ID = temp[8].replace("Parent=","")
    
    if ID in IDS:
        pass
    else:
        print i.strip()

F1.close()
