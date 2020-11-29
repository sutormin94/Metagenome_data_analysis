###############################################
##Dmitry Sutormin, 2020##
##Metagenomes analysis##

# 1) Takes metagenome assembly as a multifasta file.
# 2) Takes translated ORFs marked with MetaGeneMark and translated by InterProScan. !Tabs (\t) should be replaced beforehead by __ symbol!
# 3) Takes gff3 file with results of InterProScan annotation. !Double-semicolons (;;) and (; ) should be replaced by one semicolon and . correspondingly beforehead!
###############################################

#######
#Packages to be imported.
#######

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import Bio
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline, NcbimakeblastdbCommandline
import pandas as pd
from pandas import DataFrame
from Bio.Alphabet import generic_dna, generic_protein
from BCBio import GFF


#######
#Input data.
#######

#Path to working directory.
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Spongy\GFF_and_FASTA_to_GB\\"
#Metagenome assembly.
Metagenome_assembly_path=PWD+"INPUT_Files_required\Palmata_s3_2016_contigs_more_than_5000.fasta"
#IPS ORF translations.
ORF_translations_path=PWD+"INPUT_Files_required\H_palmata_cut_cont_2016_prot.faa"
#IPS GFF3 annotation results.
ORF_annotation_path=PWD+"INPUT_Files_required\H_palmata_cut_cont_2016_prot.faa.gff3"

#Sample definition, like: Homoeodycta palmata metagenome, bacterial fraction, 2016
Sample_definition="Homoeodycta palmata metagenome, bacterial fraction, 2016"
#Output folder.
Output_folder=PWD+"OUTPUT_BCT_H_palmata_2016_IPS_GB\\"


#######
#Reads FASTA file with reference genome
#######

def read_genome(genome_path):
    #Opens FASTA file with the metagenome assembly.
    fasta_input = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta", generic_dna))
    #for name, data in fasta_input.items():
    #    print(name, data.seq)
    return fasta_input

Contigs_dict=read_genome(Metagenome_assembly_path)

#######
#Reads FASTA file with translated ORFs. Combines ORFs belongs to one contig.
#######

def read_tORF_combine(ORF_path):
    #Opens FASTA file with ORFs translations.
    fasta_input = SeqIO.to_dict(SeqIO.parse(ORF_path, "fasta", generic_protein))   
    orfs_by_contigs={}
    for name, data in fasta_input.items():
        #Get names.
        name=name.split('__')
        ORF_name=name[0]
        Contig_name=name[1].lstrip('>')
        #Assemble ORFs by contigs.
        if Contig_name in orfs_by_contigs:
            orfs_by_contigs[Contig_name][ORF_name]=data.seq
        else:
            orfs_by_contigs[Contig_name]={ORF_name : data.seq}
        
    return orfs_by_contigs

ORFs_by_contig_dictionary=read_tORF_combine(ORF_translations_path)

#######
#Reads GFF3 file with annotation by InterProScan. Combine annotation features belongs to one ORF.
#######

def read_GFF3_combine(GFF3_path):
    GFF_file=open(GFF3_path, 'r')
    ORFs_dict={}
    start_checker=0
    for line in GFF_file:
        #Check if header or new ORF.
        if line[0] in ['#']:
            line=line.rstrip().split(' ')
            
            #Check if a new ORF strats.
            if line[0] in ['##sequence-region']:  
                #Some annotation is already ready and should be stored.
                if start_checker==1:
                    #Fill /note field.
                    info_string=''
                    if len(INFO_to_keep_dict)>0:
                        for DB, data in INFO_to_keep_dict.items():
                            for feature in data:
                                info_string=info_string+feature  
                        ORFs_dict[ORF_name]['/note']=info_string
                    
                    #Fill /product field.
                    if len(PFAM_to_keep)>0:
                        pfam_string=''
                        for feature in PFAM_to_keep:
                            pfam_string=pfam_string+feature
                        ORFs_dict[ORF_name]['/product']=pfam_string

                #Start new annotation.
                CDSs_info=line[1].split('|')
                CDS_start=CDSs_info[4]
                CDS_end=CDSs_info[5]
                CDS_strand=CDSs_info[3]
                INFO_to_keep_dict={}
                PFAM_to_keep=[]
                
                if CDS_strand=='+':
                    ORFs_dict[line[1]]={'CDS' : CDS_start+'..'+CDS_end, '/codon_start' : 1, '/experiment' : 'Predicted by MetaGeneMark and annotated by InterProScan'}
                elif CDS_strand=='-':
                    ORFs_dict[line[1]]={'CDS' : 'complement('+CDS_start+'..'+CDS_end+')', '/codon_start' : 1, '/experiment' : 'Predicted by MetaGeneMark and annotated by InterProScan'}
           
            #FASTA section of a file starts.
            elif line[0] in ['##FASTA']:
                break
        
        #Reading annotation.
        else:
            start_checker=1
            
            #Parse data.
            line=line.rstrip().split('\t')
            ORF_name=line[0]
            DB=line[1]
            feature_type=line[2]
            feature_start=line[3]
            feature_end=line[4]
            feature_info_ar=line[8].split(';')
            features_dict={}
            for feature in feature_info_ar:
                feature=feature.split('=')
                features_dict[feature[0]]=feature[1]
            
            #Convert to genebank.
            if DB in ['Pfam', 'TIGRFAM', 'ProSiteProfiles', 'PRINTS', 'SUPERFAMILY', 'Gene3D']:
                #Retrive features information.
                signame=DB + ': ' 
                if 'signature_desc' in features_dict:
                    signame=signame + features_dict['signature_desc'] + '; '
                if 'Name' in features_dict:
                    signame=signame + features_dict['Name'] + '; '
                
                #Keep info in dict.
                if DB in INFO_to_keep_dict:
                    if signame not in INFO_to_keep_dict[DB]:
                        INFO_to_keep_dict[DB].append(signame)
                else:
                    INFO_to_keep_dict[DB]=[signame]
                
                #Keep PFAM data additionally.
                if DB=='Pfam':
                    if signame not in PFAM_to_keep:
                        PFAM_to_keep.append(signame)
                          
    GFF_file.close()
    
    return ORFs_dict

ORFs_data_dictionary=read_GFF3_combine(ORF_annotation_path)


#######
#Fit strings to be shorter 80 characters.
#######

def fit_string(input_string, prefix_start):
    
    #Check if line is already fitted.
    if (len(prefix_start)+len(input_string)+len('"'))<81:
        prefix_start=prefix_start+input_string+'"\n'
    else:
        input_string_split=input_string.split(' ')
        #Split extra-long words additionally.
        input_string_split_split=[]
        for i in range(len(input_string_split)):
            if len(input_string_split[i])>60:
                #print(input_string_split[i])
                subwords_ar=input_string_split[i].split(';')
                for word in subwords_ar:
                    if len(word)>60:
                        input_string_split_split.append(word[:int(float(len(word))/2)])
                        input_string_split_split.append(word[int(float(len(word))/2):])
                    else:
                        input_string_split_split.append(word)
            else:
                input_string_split_split.append(input_string_split[i])
                
        input_string_split=input_string_split_split
        #print(input_string_split, prefix_start)
        #print(qwerty)
        
        #Complete the first line of a feature.
        while len(input_string_split)>0:
            if (len(prefix_start)+len(input_string_split[0]))<80:
                prefix_start=prefix_start+input_string_split[0]+' '
                input_string_split=input_string_split[1:]
                #print(input_string_split)
                if len(input_string_split)==0:
                    prefix_start=prefix_start+'"\n'
            else:
                break
        
        #Finish the first line.
        if len(input_string_split)>0:
            prefix_start+='\n'
        
        #print(prefix_start)
        
        
        #Transfer info onto new lines.
        prefix='                     '
        while (len(input_string_split)>0):
            #print(prefix_start)
            #print(input_string_split)
            if (len(prefix)+len(input_string_split[0]))<80:
                prefix=prefix+input_string_split[0]+' '
                input_string_split=input_string_split[1:]   
                if len(input_string_split)==0:
                    prefix_start=prefix_start+prefix+'"\n'
                else:
                    continue
            else:
                prefix_start=prefix_start+prefix+'\n'
                prefix='                     '
    
    return prefix_start


#######
#Write GeneBank file.
#######

def take_all_return_genbank(contigs_dict, orfs_dict, orfs_data_dict, Sample_def, gbk_pwd):
    
    #Iterate contigs.
    for contig_name, contig_sequence in contigs_dict.items():
        contig_length=contig_name.split('_')[3]
        #Write sequence header.
        genebank_outfile=open(gbk_pwd+contig_name+'_MGM_IPS.gb', 'w')
        
        genebank_outfile.write(f'LOCUS       {contig_name.split("_")[0]}_{contig_name.split("_")[1]}          {contig_length} bp    DNA     linear CON 13-APR-2020\n')
        genebank_outfile.write(f'DEFINITION  {Sample_def}\n')
        genebank_outfile.write('SOURCE      Environmental sample\n')
        genebank_outfile.write('FEATURES             Location/Qualifiers\n')
        genebank_outfile.write(f'     source          1..{contig_length}\n')
        genebank_outfile.write('                     /mol_type="metagenomic DNA\n')
        genebank_outfile.write('                     /country="Russia, Moscow"\n')
        genebank_outfile.write('                     /collection_date="2016"\n')
        genebank_outfile.write('                     /note="from K.V. Severinov laboratory"\n')

        #Iterate ORFs annotated in the contig.
        if contig_name in orfs_dict:
            for ORF_name, ORF_seq in orfs_dict[contig_name].items():
                
                #Get ORF info.
                if ORF_name in orfs_data_dict:
                    ORF_information=orfs_data_dict[ORF_name]
                else:
                    #If CDS was not annotated by IPS, create mock annotation.
                    CDSs_info=ORF_name.split('|')
                    CDS_start=CDSs_info[4]
                    CDS_end=CDSs_info[5]
                    CDS_strand=CDSs_info[3]   
                    if CDS_strand=='+':
                        ORF_information={'CDS' : CDS_start+'..'+CDS_end, '/codon_start' : 1, '/experiment' : 'Predicted by MetaGeneMark and annotated by InterProScan'}
                    elif CDS_strand=='-':
                        ORF_information={'CDS' : 'complement('+CDS_start+'..'+CDS_end+')', '/codon_start' : 1, '/experiment' : 'Predicted by MetaGeneMark and annotated by InterProScan'}
                #print(ORF_information)
                
                #Write CDS information.
                genebank_outfile.write(f'     CDS             {ORF_information["CDS"]}\n')
                genebank_outfile.write(f'                     /codon_start={ORF_information["/codon_start"]}\n')
                
                if "/product" in ORF_information:
                    product_string_fitted=fit_string(ORF_information["/product"], '                     /product="')
                    genebank_outfile.write(f'{product_string_fitted}')
                
                if "/note" in ORF_information:
                    note_string_fitted=fit_string(ORF_information["/note"], '                     /note="')
                    genebank_outfile.write(f'{note_string_fitted}')   
                
                experiment_string_fitted=fit_string(ORF_information["/experiment"], '                     /experiment="')
                genebank_outfile.write(f'{experiment_string_fitted}')        
               
                #Write ORF translation.
                if len(ORF_seq)<45:
                    genebank_outfile.write(f'                     /translation="{ORF_seq}"\n')
                else:
                    genebank_outfile.write(f'                     /translation="{ORF_seq[:44]}\n')
                    seq_remainder=ORF_seq[44:]
                    num_of_strings=len(seq_remainder)//58
                    if num_of_strings==0:
                        genebank_outfile.write(f'                     {seq_remainder}"\n')
                    else:
                        for i in range(num_of_strings):
                            genebank_outfile.write(f'                     {seq_remainder[i*58:((i+1)*58)]}\n')
                        genebank_outfile.write(f'                     {seq_remainder[num_of_strings*58:]}"\n')
        
        #Write ORIGIN sequence - contig sequence.
        genebank_outfile.write('ORIGIN      \n')
        
        contig_sequence=str(contig_sequence.seq)
        num_of_strings=len(contig_sequence)//60
        for i in range(num_of_strings):
            start=(i*60)+1
            num_blanks=9-len(str(start))
            genebank_outfile.write(f'{" "*num_blanks}{start}')
            for j in range(6):
                genebank_outfile.write(f' {contig_sequence[((start-1)+(j*10)):((start-1)+((j+1)*10))]}')
            genebank_outfile.write('\n')
        
        #Write remainder of a sequence.
        start=start+60
        num_blanks=9-len(str(start))
        genebank_outfile.write(f'{" "*num_blanks}{start}')
        
        seq_remainder=contig_sequence[start-1:]
        num_of_blocks_rem=len(seq_remainder)//10
        for k in range(num_of_blocks_rem):
            start_rem=k*10
            genebank_outfile.write(f' {seq_remainder[start_rem:(start_rem+10)]}')
        genebank_outfile.write(f' {seq_remainder[num_of_blocks_rem*10:]}\n//\n')
           
        genebank_outfile.close()
    
    return

take_all_return_genbank(Contigs_dict, ORFs_by_contig_dictionary, ORFs_data_dictionary, Sample_definition, Output_folder)