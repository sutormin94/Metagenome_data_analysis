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
Taxonomy_MGRAST="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\MGRAST_strains_taxonomy\\NTC.csv"

#NCBI taxonomy names (convert synonims).
NCBI_tax_names="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\\NCBI_taxonomy\\names.dmp"

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
        #Count file lines.
        i=0
        #Keep species names with abundance.
        Species_dict={}
        for line in filein:
            line=line.rstrip().split('\t')
            #Track progress.
            i+=1
            if i%10000==0:
                print(str(i)+' lines processed.')
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
        Species_dict_sorted=Species_dict
        Species_dict_sorted={k: v for k, v in sorted(Species_dict_sorted.items(), key=lambda item: item[1], reverse=True)}        
        print('Number of species identified in metagenome: ' + str(len(Species_dict)))
    
    elif db_format=='MGRAST':
        filein=open(fileinpath, 'r')
        #Count file lines.
        i=0  
        #Keep species names with abundance.
        Species_dict={}
        Strains_dict={}        
        for line in filein:
            line=line.rstrip().split('\t')
            #Track progress.
            i+=1
            if i%10000==0:
                print(str(i)+' lines processed.')            
            #Filter only species-level-resolved.
            if line[-1] != "semicolon separated list of annotations":
                spec_data=line[-1].split(';')
                #print(spec_data)
                for strain_id in spec_data:
                    strain_id=strain_id.rstrip(']').lstrip('[')
                    #Keep all strains with frequences.
                    if strain_id in Strains_dict:
                        Strains_dict[strain_id]+=1
                    else:
                        Strains_dict[strain_id]=1
                    
                    spec_id=strain_id.split(' ')
                    #Keep all species with frequences.
                    #print(spec_id)
                    if len(spec_id)>1:
                        if spec_id[1] != 'sp.':
                            spec_id=spec_id[0]+' '+spec_id[1]
                            if spec_id in Species_dict:
                                Species_dict[spec_id]+=1
                            else:
                                Species_dict[spec_id]=1

        filein.close()  
        Species_dict_sorted=Species_dict
        Species_dict_sorted={k: v for k, v in sorted(Species_dict_sorted.items(), key=lambda item: item[1], reverse=True)}
        print('Number of species identified in metagenome: ' + str(len(Species_dict)))
        print('Number of strains identified in metagenome: ' + str(len(Strains_dict)) + '\n\n')
    return Species_dict, Species_dict_sorted

#Dictionary of species identified in metagenomic data with its abundance.
Spec_dict,Species_dict_sorted=read_met_taxo(Taxonomy_MGRAST, 'MGRAST')
#Write the data.
Write_taxo_data(Spec_dict, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\MGRAST_strains_taxonomy\\NTC_MGRAST_metagenome_taxonomy_species_only')
Write_taxo_data(Species_dict_sorted, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\MGRAST_strains_taxonomy\\NTC_MGRAST_metagenome_taxonomy_species_only_sorted')


##############
########  Part 2. Working with Ezbiocloud taxonomy, NCBI taxons names to detect synonimous names, 
########  retriving species names and intersecting with metagenomes data obtained in Part 1.
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
#Read NCBI taxon names, make dictionary - Name : NCBI_tax_id.
#######

def read_ncbi_tax_names(NCBI_tax_names_path):
    filein=open(NCBI_tax_names_path, 'r')
    name_id_dict={}
    id_syn_names_dict={}
    count_amb=0
    i=0
    for line in filein:
        line=line.rstrip().split('\t')
        #Track progress.
        i+=1
        if i%10000==0:
            print(line[0], line[2], line)
            print(str(i)+' lines processed.')         
        #Construct name_id_dict.
        if line[2] not in ['environmental samples']:
            if line[2] in name_id_dict:
                print('Already here: ambigous name - ' + line[2])
                count_amb+=1            
            name_id_dict[line[2]]=line[0]
            
        
        #Construct id_syn_names_dict.
        if line[0] in id_syn_names_dict:
            id_syn_names_dict[line[0]].append(line[2])
        elif line[0] not in id_syn_names_dict:
            id_syn_names_dict[line[0]]=[line[2]]
            
    filein.close()
    print('Number of taxo names detected in NCBI db: ' + str(len(name_id_dict)))
    print('Number of ambigous taxo names detected in NCBI db: ' + str(count_amb))
    print('Number of taxo_ids detected in NCBI db: ' + str(len(id_syn_names_dict)) + '\n\n')
    return name_id_dict, id_syn_names_dict

NCBI_name_id_dict, NCBI_id_syn_names_dict=read_ncbi_tax_names(NCBI_tax_names)

#######
#Find Ezbiocloud taxonomic id of metagenome species.
#######

def return_ezbio_id(Spec_dict, Ezbio_species_dict, name_id_dict, id_syn_names_dict, fileoutpath):
    #Keep identified in Ezbiocloud DB species.
    count_ident_self=0
    count_ident_by_syn=0
    Dict_of_identified={}
    #Keep not identified species.
    fileout=open(fileoutpath, 'w')
    i=0
    for species, abund in Spec_dict.items():
        #Track progress.
        i+=1
        if i%100==0:
            print(str(float(i)*100/len(Spec_dict))+'% species processed.')         
        if species in Ezbio_species_dict:
            #Identify species in Ezbiocloud db.
            taxo_id=Ezbio_species_dict[species]
            Dict_of_identified[taxo_id]=[species, abund]
            count_ident_self+=1
        else:
            if species in name_id_dict:
                #Use NCBI taxonomy names db to try synonimous names.
                species_tax_id=name_id_dict[species]
                species_syn_list=id_syn_names_dict[species_tax_id]
                ambiguity=0
                for species_syn in species_syn_list:
                    if species_syn in Ezbio_species_dict:
                        taxo_id=Ezbio_species_dict[species_syn]
                        Dict_of_identified[taxo_id]=[species_syn, abund]
                        count_ident_by_syn+=1     
                        ambiguity+=1
                if ambiguity>1:
                    print('Species name ' + species + ' is ambiguos: counted ' + str(ambiguity) + ' times.')
                elif ambiguity==0:
                    print(species, species_tax_id, species_syn_list)
                    print('Species name ' + species + ' is not matched with Ezbiocloud db even by using NCBI name synonims.')
                    fileout.write(str(abund) + '\t' + species + '\n')
            else:
                print('Species name ' + species + ' is not found in NCBI name synonims db.')
    fileout.close()
    print(count_ident_self, count_ident_by_syn, len(Spec_dict))
    print('Fraction of metagenomic species correctly identified in Ezbiocloud: ' + str(float(count_ident_self+count_ident_by_syn)/len(Spec_dict)))
    print('Species not identified are likely to have synonimus names \n\n')
    return Dict_of_identified

Dict_of_identified_in_Ezbio=return_ezbio_id(Spec_dict, Ezbio_species_dict, NCBI_name_id_dict, NCBI_id_syn_names_dict, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\MGRAST_strains_taxonomy\\NTC_MGRAST_metagenome_taxonomy_species_only_not_identified_in_Ezbiocloud')


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

return_v34_for_met_spec(Ezbio_v34_dict, Dict_of_identified_in_Ezbio, 'C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\MGRAST_strains_taxonomy\\NTC_MGRAST_metagenome_taxonomy_species_only_V3V4_seqs_from_Ezbiocloud.fasta')

############
###### Now sequences are ready to be clusterized with MMSeqs2.
############

