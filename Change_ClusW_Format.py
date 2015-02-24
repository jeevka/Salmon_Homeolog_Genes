import sys
import re

########################################################
# Pick Full Seq names
########################################################

#Seq3 = Seq1.replace("S","salmon")
#Seq4 = Seq2.replace("SH","salmonH")

########################################################
# Add Full Seq names to Alignment file
########################################################
ALN = open("Homeo_2015_18.fasta.aln","r")

#S1 = Seq1.split("-")[0]
#S2 = Seq2.split("-")[0]

for i in ALN:
    if re.search("^S-",i):
        temp = i.split()
        Seq1 = temp[0].replace("S","salmon")
        Seq1 = Seq1.replace("_+_","(+)")
        #print Seq1
        print Seq1,"\t",temp[1]
    elif re.search("^SH-",i):
        temp = i.split()
        Seq2 = temp[0].replace("SH","salmonH")
        Seq2 = Seq2.replace("_+_","(+)")
        #print Seq2
        print Seq2,"\t\t",temp[1]
    else:
        print "\t\t\t\t\t",i.strip()
        
ALN.close()