###############################################
##Dmitry Sutormin, 2019##
##16S data analysis##

#1) Script parses taxonomy for metagenome reads assigned by Ghost Koala or Kaiju.
#2) Keeps only taxonomy resolved till species level. 
#3) Retuns 16S sequences for specieses.
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

#Metagenome taxonomy data.
Taxonomy_kaiju="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\NTC_Kaiju_metagenome_taxonomy_kaiju_taxonpaths.txt"
Taxonomy_ghost_koala="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Snow_Ghost_koala_metagenome_taxonomy_user_out_topTaxonomy"
Taxonomy_KEGG="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\NTC_KEGG_metagenome_taxonomy_kaiju_taxonpaths_user_out.top"

#Ezbiocload taxonomy.
Ezbio_taxonomy="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Ezbiocloud_database\ezbiocloud_id_taxonomy.txt"

#Ezbiocloud V3-V4 sequences.
Ezbio_seqs="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Ezbiocloud_database\\V3_V4_region_from_ezbiocloud_alignment.fasta"

##############
########  Part 1. Working with metagenome taxonomy, retriving species-resolved data.
##############

#######
#Write taxonomic data.
#######

def Write_taxo_data(Species_dict, fileoutpath):
    fileout=open(fileoutpath, 'w')
    for taxo, abund in Species_dict.items():
        fileout.write(str(abund) + '\t' + taxo + '\n')
    fileout.close()
    return


#######
#Read taxonomy data from metagenome, filter species-resolved only.
#######

def read_met_taxo(fileinpath, db_format):
    if db_format=='kaiju':
        filein=open(fileinpath, 'r')
        count_long=0
        count_sp1=0
        count_spec=0
        count_complex=0
        #Keep species names with abundance.
        Species_dict={}
        for line in filein:
            line=line.rstrip().split('\t')
            #Filter only species-level-resolved.
            if len(line)>9:
                line_end_split=line[-1].split(' ')
                line_semiend_split=line[-2].split(' ')
                species_name=''
                if len(line_end_split)>2:
                    if len(line_semiend_split)>2:
                        count_long+=1
                        #print('Both names are very long!  ' + str(count_long) + '  ' + str(line[-3:]))
                        if line[-1].find('subsp.')>-1:
                            species_name=line_end_split[0] + ' ' + line_end_split[1]
                            #print('Subspecies identified!  ' + str(count_long) + '  ' + str(species_name))
                        elif line[-1].find('bacterium')>-1:
                            #print('Likely to be smth sp.!  ' + str(count_long) + '  ' + str(line[-1]))
                            continue
                    elif len(line_semiend_split)==2:
                        species_name=line[-2]
                    else:
                        count_sp1+=1
                        if line[-1].find('sp.')>-1:
                            #print('Likely to be smth sp. 1!  ' + str(count_sp1) + '  ' + str(line[-1]))
                            continue
                        else:
                            if (line[-1].find('complex')>-1) or (line[-1].find('group')>-1) or (line[-1].find('al.')>-1):
                                species_name=line_end_split[0] + ' ' + line_end_split[1]
                            elif (line[-1].find('Candidatus')>-1):
                                species_name=line_end_split[1] + ' ' + line_end_split[2]
                            else:
                                #print('What is this 1?  ' + str(count_sp1) + '  ' + str(line[-1]))
                                continue
                        
                elif len(line_end_split)==2:
                    if len(line_semiend_split)>2:
                        count_complex+=1
                        species_name=line[-1]
                        #print('Likely to be complex indication!  ' + str(count_complex) + '  ' + str(line[-3:])) 
                    elif len(line_semiend_split)==2:
                        line[-1]=re.sub('[\[\]]', '', line[-1])
                        species_name=line[-1]
                        #print('What is this?  ' + '  ' + str(line[-1])) 
                    else:
                        count_spec+=1
                        species_name=line[-1]
                        #print('Likely to be a species!  ' + str(count_spec) + '  ' + str(line[-1]))   
                #Keep final species names.
                if (species_name!='') and (species_name!='environmental samples') and (species_name.find('unclassified')==-1):
                    species_name=re.sub('[\[\]]', '', species_name)
                    if species_name in Species_dict:
                        Species_dict[species_name]=Species_dict[species_name]+int(line[0])
                    else:
                        Species_dict[species_name]=int(line[0]) 
                
        filein.close()
        print('Number of species identified in metagenome: ' + str(len(Species_dict)))
    return Species_dict

#Dictionary of species identified in metagenomic data with its abundance.
Spec_dict=read_met_taxo(Taxonomy_kaiju, 'kaiju')
#Write the data.
Write_taxo_data(Spec_dict, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\NTC_Kaiju_metagenome_taxonomy_kaiju_taxonpaths_species_only_end')


##############
########  Part 2. Working with Ezbiocloud taxonomy, retriving species names and intersecting with metagenomes data obtained in Part 1.
##############


#######
#Read taxonomy data from Ezbiocloud.
#######

def read_ezbiocloud(dbinfilepath):
    filein=open(dbinfilepath, 'r')
    #Keep data in dictionary.
    DB_dict={}
    for line in filein:
        species=line.rstrip().split('\t')[1].split(';')[-1]
        DB_dict[species]=line.rstrip().split('\t')[0]
    filein.close()
    return DB_dict

Ezbio_species_dict=read_ezbiocloud(Ezbio_taxonomy)


#######
#Find Ezbiocloud taxonomic id of metagenome species.
#######

def return_ezbio_id(Spec_dict, Ezbio_species_dict, fileoutpath):
    #Keep identified in Ezbiocloud DB species.
    count_ident=0
    Dict_of_identified={}
    #Keep not identified species.
    fileout=open(fileoutpath, 'w')
    for species, abund in Spec_dict.items():
        if species in Ezbio_species_dict:
            taxo_id=Ezbio_species_dict[species]
            Dict_of_identified[taxo_id]=[species, abund]
            count_ident+=1
        else:
            fileout.write(str(abund) + '\t' + species + '\n')
    fileout.close()
    print('Fraction of metagenomic species correctly identified in Ezbiocloud: ' + str(float(count_ident)/len(Spec_dict)))
    print('Species not identified are likely to have synonimus names')
    return Dict_of_identified

Dict_of_identified_in_Ezbio=return_ezbio_id(Spec_dict, Ezbio_species_dict, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\NTC_Kaiju_metagenome_taxonomy_kaiju_taxonpaths_species_only_not_identified_in_Ezbiocloud')


##############
########  Part 3. Returning V3-V4 sequences of metagenome species matched with Ezbio DB in Part 2.
##############


#######
#Reads mfa with V3-V4, returns dict of sequences.
#######

def read_mfa_make_dict(fileinpath):
    filein=open(fileinpath, 'r')
    #Keep seq data in dict.
    V34_dict={}
    i=0
    for record in SeqIO.parse(filein, "fasta"):
        seq_v3v4=str(record.seq)
        id_v3v4=record.name
        V34_dict[id_v3v4]=seq_v3v4
        i+=1
    filein.close()  
    print('Number of V3-V4 sequences avaliable in Ezbio DB: ' + str(len(V34_dict)))
    return V34_dict

Ezbio_v34_dict=read_mfa_make_dict(Ezbio_seqs)


#######
#Return V3-V4 sequences for species from metagenome matched with Ezbio DB.
#######

def return_v34_for_met_spec(V34_dict, met_spec_dict, fileoutpath):
    fileout=open(fileoutpath, 'w')
    count_v34=0
    for tax_id, data in met_spec_dict.items():
        if tax_id in V34_dict:
            #Fasta names: TaxID_Species name_Abundance
            fileout.write('>'+tax_id+'_'+data[0]+'_'+str(data[1])+'\n'+V34_dict[tax_id]+'\n')
            count_v34+=1
        else:
            print('Smth wrong: no sequence with id provided - ' + tax_id)
    fileout.close()
    print('Number of metagenome species with V3-V4 sequences avaliable in Ezbio DB: ' + str(count_v34))
    return

return_v34_for_met_spec(Ezbio_v34_dict, Dict_of_identified_in_Ezbio, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\NTC_Kaiju_metagenome_taxonomy_kaiju_taxonpaths_species_only_V3V4_seqs_from_Ezbiocloud.fasta')

############
###### Now sequences are ready to be clusterized with MMSeqs2.
############

