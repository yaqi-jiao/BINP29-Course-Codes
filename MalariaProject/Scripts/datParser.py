
# -*- coding: utf-8 -*-

"""
This script is created to retrieve scaffolds that may contain host genes, by comparing the BLAST-hits with the swissprot database.

Runs like:

python datParser.py Ht.blastx Ht.fna taxonomy.dat uniprot_sprot.dat

*********

Approach:

* Retrieve the first hit for each query from the blast file, and put these names in a dictionary together with the query ID.

* Parse through uniprot dat-file to retrieve GeneID number to each corresponding accessionnumber (hits)

* Look in the taxfile if the geneID is found in there, and store the queries in a list/set.

* Go through the query list and parse through the gff parsed fasta file to retrieve which scaffolds these correspond to.

* Print scaffolds that derive from bird.
"""

#%%

import sys

try:
    blastfile = sys.argv[1]
    fastafile = sys.argv[2]
    taxfile = sys.argv[3]
    uniprot = sys.argv[4]
except:
    print("Run program like: datParser.py Ht.blastx Ht.fna taxonomy.dat uniprot_sprot.dat")
    sys.exit()
#%%
################## PARSE BLAST FILE #############################
# The Q boolean is to control that only the first hit for each query is included.
Q = False
AC_dict = {}
with open(blastfile, "r") as BLAST:
    for line in BLAST:
        if line.startswith("Query= "):
            Query = line.rstrip().split()[-1]
            Q = True
        elif line.startswith(">"):
            # If Q is true, means that we have found a Query previously.
            if Q == True:
                # Different blast versions outputs different Hit format. Some like >HIT_ID.x some like > sp|HIT_ID|, and some like >sp|HIT_ID|.
                # Following will check if | in line, and split with respect to that, if not, it will expect format like >HIT_ID.x
                if "|" in line:
                    AC = line.split("|")[1]
                else:
                    line = line.split()[0]
                    AC = line.split(".")[0][1:]

                # Now add Query and AC to dict.
                AC_dict[AC] = [Query]
                # Set Q to false, so that if the query has multiple hits, we wont add these.
                Q = False

###################################################################
#%%

##################### Get GeneIDset for birdtaxonomy ##############
aves_found = False
SN = ''
SNset = set()
with open(taxfile, "r") as T1:
    for line in T1:
        # Ensuring we are in the "birds" section of the taxonomy file.
        if aves_found == True:
            # We want to keep the name from "SCIENTIFIC NAME" row
            if line.startswith("SCIENTIFIC"):
                if SN:
                    SNset.add(SN)
                SN = line.rstrip().split(":")[-1][1:]

            elif line.startswith("BLAST"):
                aves_found = False
                break
        # This part here will set the aves_found to True, when the BLAST NAME "birds" is found in the taxonomy file.
        # All rows below will be species belonging to bird, until next "BLAST NAME" is found.
        elif line.startswith("BLAST NAME"):
            line = line.split(":")[1]
            if "birds" in line:
                aves_found = True


#%%

def check(SNlist):
    # Function to check if any of the items in the list are also found in the SNset - meaning that the item is a bird species.
    for S in SNlist:
        if S in SNset:
            return(True)

    return(False)

##################### Parse uniprot_sprot.dat #####################
birdIDs = set()
OC = ''
AC = ''
i=0
OC_dict = {}
AClist = []
with open(uniprot, "r") as UNIPROT:
    # Parse through the uniprot_sprot.dat file, looking at the AC row to identify accessionnumbers from the Blast query dictionary.
    # for each AC, retrieve the rows beginning with OC - which contains species information that also exist in the taxonomy.dat file
    for line in UNIPROT:
        if line.startswith("AC"):
            if AC:
                if AC in AC_dict:
                    if check(OC) == True:
                        AClist += AC_dict[AC]
            AC = line.rstrip().split()[-1][:-1]
            OC = ''
        elif line.startswith("OC"):
            line = line.rstrip().replace(" ", "")[2:-1]
            line = line.split(";")
            if OC:
                OC += line
            else:
                OC = line


##################################################################
#%%

#################### Retrieve Scaffold list ######################
                
                
with open(fastafile, "r") as Fasta:
    for line in Fasta:
        if line.startswith(">"):
            ID = line.split()[0][1:]
            if ID in AClist:
                print(line.split()[2].split("=")[-1])
    
###################################################################
    
    
