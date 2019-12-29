###############################################
##Dmitry Sutormin, 2019##
##16S data analysis##

#Script subtracts a fragment of 16S rRNA gene from multiple alignment by coordinates provided (V3-V4 region by default).
###############################################


#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
import re
import matplotlib.pyplot as plt


#######
#Data to be used.
#######

#Input 16S genes multiple alignment (EzBioCloud release 12.2019).
#Alig_16S="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Ezbiocloud_database\ezbiocloud_full_align.fasta"
Alig_16S="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Grotta_sponges_seqtab_nochim_mergers.fasta"

#Output file with sequences.
Alig_16S_subtr="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Grotta_sponges_seqtab_nochim_mergers_primer_trimmed.fasta"

#Coordinates of subtraction.
#Left_start=9238-1
#Left_end=9282
#Right_start=14969-1
#Right_end=15349
Left_end=0
Right_start=-4

#######
#Read multiple alignment subtract fragment.
#######

def read_subtract(fileinpath, fileoutpath, coord_l, coord_r):
    seq_len_ar=[]
    filein=open(fileinpath, 'r')
    fileout=open(fileoutpath, 'w')
    i=0
    for record in SeqIO.parse(filein, "fasta"):
        if i%10000 == 0:
            print(str(i) + ' done')
        seq_16S=str(record.seq)[coord_l:coord_r]
        seq_16S_red=re.sub('[-]', '', seq_16S)
        id_16S=record.name
        fileout.write('>'+id_16S+'\n'+seq_16S_red+'\n')
        seq_len_ar.append(len(seq_16S_red))
        i+=1
    filein.close()  
    fileout.close()
    return seq_len_ar

Seq_len_array=read_subtract(Alig_16S, Alig_16S_subtr, Left_end, Right_start)

#Plot distribution of the length of V3-V4 region.
plt.hist(Seq_len_array)
plt.yscale('log')
plt.xlabel('Length of V3-V4 regions from Grotta merged dataset')
plt.ylabel('Number of sequences')
plt.show()
plt.savefig('C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\Grotta_sponges_seqtab_nochim_mergers_primer_trimmed.png')