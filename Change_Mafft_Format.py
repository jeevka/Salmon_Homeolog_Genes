import sys
import re

########################################################
# Pick Full Seq names
########################################################
Seq = open("Homeo_2015_11.fasta","r")

N = 0
for i in Seq:
    if re.search(">",i):
        N += 1
        if N == 1:
            Seq1 = (i.split(">")[1]).strip()
        else:
            Seq2 = (i.split(">")[1]).strip()

Seq.close()


Seq3 = Seq1.replace("S","salmon")
Seq4 = Seq2.replace("SH","salmonH")

########################################################
# Add Full Seq names to Alignment file
########################################################
ALN = open("/mnt/users/jeevka/Salmon_Homeolog_Alignments/Mafft_Algnment_Results/Homeo_2015_11.aln","r")

S1 = Seq1.split("-")[0]
S2 = Seq2.split("-")[0]

for i in ALN:
    if re.search("^S-",i):
        temp = i.split()
        print Seq3,"\t",temp[1]
    elif re.search("^SH-",i):
        temp = i.split()
        print Seq4,"\t",temp[1]
    else:
        print "\t\t\t\t\t",i.strip()
        
ALN.close()
