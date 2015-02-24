import sys
import re

# Change the range based on strand
def check_for_reverse_strand_1(t1,t2):
    if int(t1) > int(t2):
        H_R1 = int(t2)
        H_R2 = int(t1)
    else:
        H_R1 = int(t1)
        H_R2 = int(t2)
        
    return H_R1,H_R2

# Change the range based on strand
def check_for_reverse_strand_2(t1,t2,strand,N):
    if int(t1) > int(t2):
        H_R1 = int(t2)
        H_R2 = int(t1)
        strand[N] = "-"
    else:
        H_R1 = int(t1)
        H_R2 = int(t2)
        strand[N] = "+"
        
    return H_R1,H_R2,strand

def find_new_positions_1(OR1,OR2,HR1,HR2,R1,R2):
    # Distance from Homeolog start
    HD1 = int(R1) - int(HR1) + 1 
    HD2 = int(R2) - int(HR1) + 1
    
    # New starting point for Homeolog Genes
    new_start = OR1 + HD1
    new_end = OR1 + HD2
    
    return new_start,new_end

def find_new_positions_2(OR1,OR2,HR1,HR2,R1,R2):
    # Distance from Homeolog start
    HD1 = int(HR1) - int(R1) + 1 
    HD2 = int(HR1) - int(R2) + 1
    
    # New starting point for Homeolog Genes
    new_start = OR1 + HD1
    new_end = OR1 + HD2
    
    return new_start,new_end

def get_function_information(temp):
    note ="" #temp[8].split("=")[-1]
    for i in xrange(9,len(temp)):
        note += " " + temp[i]
    
    return note

def edit_feature_1(temp):
    temp1 = temp.split(";")
    text = ""
    N = 0
    for i in temp1:
        N += 1
        if re.search("^ID",i) or re.search("^Name",i):
            i = i + "_Homeolog"
        if N == len(temp1):
            text = text + i
        else:
            text = text + i + ";"
        
    return text                

def edit_feature_2(temp):
    temp1 = temp.split(";")
    text = ""
    N = 0
    for i in temp1:
        N += 1
        if re.search("^Parent",i):
            i = i + "_Homeolog"
        if N == len(temp1):
            text = text + i
        else:
            text = text + i + ";"
        
    return text
    
def pick_homeolog_genes_1(GFF3,OChr,HChr,OR1,OR2,HR1,HR2,strand,N):
    mRNA = 0
    #print OChr,OR1,OR2
    for i in GFF3:
        temp = i.split()
        # Check whether the gene models falls within the Homeolog regions
        if temp[0] == HChr and int(temp[3]) >= int(HR1) and int(temp[4]) <= int(HR2):
            if temp[2] == "mRNA" or temp[2] == "gene" or temp[2] == "CDS" or temp[2] == "five_prime_UTR" or temp[2] == "three_prime_UTR":

                # notes are multiple words
                # Get the funtion write from Notes
                if len(temp) > 9:
                    note = get_function_information(temp)
                else:
                    note = ""
                
                alias = HChr + ":" + temp[3] + "-" + temp[4]
                    
                # Calculate the new co-ordinates
                if strand[N] == "+":
                    new_start,new_end = find_new_positions_1(OR1,OR2,HR1,HR2,temp[3],temp[4])
                else:
                    new_start,new_end = find_new_positions_2(OR1,OR2,HR1,HR2,temp[3],temp[4])
                    
                if temp[2] == "gene":
                    txt = OChr + "\t" + temp[1] + "\t" + "Homeolog_Genes" + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + temp[8] + note +";Derives_from=" + alias
                
                if temp[2] == "mRNA":    
                    # 9th coloumn in gff3 file needs more changes
                    feature_text = edit_feature_1(temp[8])
                    mRNA = 1
                    txt = OChr + "\t" + temp[1] + "\t" + "Homeolog_Transcript" + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + feature_text + note + ";Derives_from=" + alias
                    
                if (temp[2] == "five_prime_UTR" or temp[2] == "three_prime_UTR" or temp[2] == "CDS") and mRNA == 1:
                    # 9th coloumn in gff3 file needs more changes
                    feature_text = edit_feature_2(temp[8])                    
                    txt = OChr + "\t" + temp[1] + "\t" + temp[2] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + feature_text
                
                if temp[2] == "sequence_assembly" or temp[2] == "three_prime_UTR":
                    mRNA = 0
                
                # In case if the CDS does not have mRNA
                # We cant load CDS features in GFF3 without mRNA (parent).
                try:
                    print txt
                except:
                    pass
            
    return 0

def pick_homeolog_genes_2(GFF3,OChr,HChr,OR1,OR2,HR1,HR2,strand,N):
    mRNA = 0
    for i in GFF3:
        temp = i.split()
         
        # Check whether the gene models falls within the Homeolog regions
        if temp[0] == HChr and int(temp[3]) >= HR1 and int(temp[4]) <= HR2:
            if temp[2] == "mRNA" or temp[2] == "gene" or temp[2] == "CDS" or temp[2] == "five_prime_UTR" or temp[2] == "three_prime_UTR":
                
                # notes are multiple words
                # Get the funtion write from Notes
                if len(temp) > 9:
                    note = get_function_information(temp)
                else:
                    note = ""
                
                alias = HChr + ":" + temp[3] + "-" + temp[4]
                    
                # Calculate the new co-ordinates
                if strand[N] == "+":
                    new_start,new_end = find_new_positions_1(OR1,OR2,HR1,HR2,temp[3],temp[4])
                else:
                    new_start,new_end = find_new_positions_2(OR1,OR2,HR1,HR2,temp[3],temp[4])

                if temp[2] == "gene":
                    txt = OChr + "\t" + temp[1] + "\t" + "Homeolog_Genes" + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + temp[8] + note +";Derives_from=" + alias
                
                if temp[2] == "mRNA":
                    # 9th coloumn in gff3 file needs more changes
                    feature_text = edit_feature_1(temp[8])
                    mRNA = 1
                    txt = OChr + "\t" + temp[1] + "\t" + "Homeolog_Transcript" + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + feature_text + note + ";Derives_from=" + alias
                    
                if (temp[2] == "five_prime_UTR" or temp[2] == "three_prime_UTR" or temp[2] == "CDS") and mRNA == 1:
                    # 9th coloumn in gff3 file needs more changes
                    feature_text = edit_feature_2(temp[8])                    
                    txt = OChr + "\t" + temp[1] + "\t" + temp[2] + "\t" + str(new_start) + "\t" + str(new_end) + "\t" + temp[5] + "\t" + temp[6] + "\t" + temp[7] + "\t" + feature_text
                    
                if temp[2] == "sequence_assembly" or temp[2] == "three_prime_UTR":
                    mRNA = 0
                
                # In case if the CDS does not have mRNA
                # We cant load CDS features in GFF3 without mRNA (parent).
                try:
                    print txt
                except:
                    pass
            
    return 0
############################################################################################
################################## MAIN PROGRAAM ###########################################
############################################################################################

# Latest GFF3 File
GFF3 = open("/mnt/users/jeevka/Make_GFF3_For_GBrowse/Salmon_3p6_Chr_051214.gff3","r")
#GFF3 = GFF3_1.readlines()

# Homeolog File
F1 = open("Ssa-CIGENE36-homologeous-regions.txt","r")

strand = {}
N = 0
for i in F1:
    temp = i.split()
    N += 1
    # Original Chromosome
    OChr = temp[0]
    
    # Homeolog Chromsome
    HChr = temp[3]
    
    # Ranges from Original region
    OR1,OR2= check_for_reverse_strand_1(temp[1],temp[2])
    
    # Ranges from "Homeolog" region
    HR1,HR2,strand = check_for_reverse_strand_2(temp[4],temp[5],strand,N)

    # Input data contains      
    # Pick the Gene models in the Homeolog regions
    pick_homeolog_genes_1(GFF3,OChr,HChr,OR1,OR2,HR1,HR2,strand,N)
    GFF3.seek(0)
    
    # Pick the Gene models in the Homeolog regions
    pick_homeolog_genes_2(GFF3,HChr,OChr,HR1,HR2,OR1,OR2,strand,N)
    GFF3.seek(0)
    
F1.close()
