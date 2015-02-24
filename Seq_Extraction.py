import sys
import re
import os
import subprocess

#################################################################################################################
####################################### SUB PROGRAM ############################################################
#################################################################################################################
def split_ranges(HR1,HR2,HChr):
    ranges = []
    if int(HR1) < int(HR2):
        R1 = HR1
        R2 = HR2
    else:
        R1 = HR2
        R2 = HR1
        
    N = int(R1)
    while N <= int(R2):
        N += 500000
        Seq = HChr+ ":" + str(R1) + "-" + str(N)
        R1 = N + 1 
        print Seq
    
    sys.exit()

def find_strand(HChr,HR1,HR2):
    # Extracting Homeolog Sequences
    if int(HR1) > int(HR2):
        Seq = HChr+ ":" + str(HR2) + "-" + str(HR1)
        strand = "-"
    else:
        Seq = HChr+ ":" + str(HR1) + "-" + str(HR2)
        strand = "+"
        
    return Seq, strand

def change_format_1(File1,F2,text,Chr,strand):
    F1 = open(File1,"r")
    F2 = open(F2,"w")
    for i in F1:
        if re.search("^>",i):
            temp = i.split(":")
            txt = ">" + text + "-" + Chr + "(" + strand + ")/" + temp[1]
            F2.write(txt)
        else:
            F2.write(i)
    
    F1.close()
    F2.close()

def change_format_2(File1,F2,text,Chr,strand):
    F1 = open(File1,"r")
    F2 = open(F2,"w")
    for i in F1:
        if re.search("^>",i):
            temp = i.split(">")
            txt = ">" + text + "-" + Chr + "(" + strand + ")/" + temp[1]
            F2.write(txt)
        else:
            F2.write(i)
    
    F1.close()
    F2.close()
#################################################################################################################
####################################### MAIN PROGRAM ############################################################
#################################################################################################################
File1 = sys.argv[1]
File2 = sys.argv[2]
File3 = sys.argv[3]

#os.system("module load samtools")
#os.system("module load clustalw/2.1")

F1 = open("Salmon_Homeolog_Regions_130215.csv","r")
N = 0
for i in F1:
    N += 1
    Fname = "Homeo_2015_" + str(N) + ".fasta"
    temp = i.split()
    OChr = temp[0]; R1 = temp[1]; R2 = temp[2]
    HChr = temp[3]; HR1 = temp[4]; HR2 = temp[5]
    
    # Extracting Sequences
    # Extracting Original Sequencese
    # Seq = OChr+ ":" + str(R1) + "-" + str(R2)
    
    Seq, strand = find_strand(OChr,R1,R2)
    
    CMD1 = "samtools faidx /mnt/users/jeevka/Salmon_3p6Chr/Salmon_3p6Chr_Unmasked.fasta " + Seq + " >" + File1
    os.system(CMD1)
    
    # Change format for Synteny browser
    change_format_1(File1,"File_1_1.fasta","S",OChr,strand) 
    
    # Extracting Homeolog Sequence
    Seq, strand = find_strand(HChr,HR1,HR2)
    
    # In the latest file, strand information was given
    strand = temp[6]
    
    CMD2 = "samtools faidx /mnt/users/jeevka/Salmon_3p6Chr/Salmon_3p6Chr_Unmasked.fasta " + Seq + " >" + File2
    os.system(CMD2)    
    
    if strand == "-":    
        cmd_2 = "revseq -reverse Y -complement N -sequence " + File2 +  " -outseq Seq_Rev.fasta"
        os.system(cmd_2)
        print "Reverse"
        # Change format for Synteny browser
        change_format_2("Seq_Rev.fasta","File_2_1.fasta","SH",HChr,strand)
        CMD3 = "cat " + "File_1_1.fasta" + " " +  "File_2_1.fasta" + " >" + Fname
    else:
        # Change format for Synteny browser
        change_format_1(File2,"File_2_1.fasta","SH",HChr,strand)        
        CMD3 = "cat " + "File_1_1.fasta" + " " +  "File_2_1.fasta" + " >" + Fname
    
    os.system(CMD3)