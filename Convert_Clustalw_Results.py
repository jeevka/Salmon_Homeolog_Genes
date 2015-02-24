import sys
import re

F1 = open("Test_1000.fasta.aln","r")

for i in F1:
    if re.search("^ssa",i):
        temp1 = i.split()
        temp2 = temp1[0].split("_")
        txt = "salmon_homeolog-"+ temp2[0] + "(" + "+" + ")/" + temp2[1] + "\t" + temp1[1]
        print txt
        sys.exit()
    else:
        print i.strip()