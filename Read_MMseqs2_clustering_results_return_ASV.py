###############################################
##Dmitry Sutormin, 2019##
##16S data analysis##

#1) Script parses MMseqs2 clustering results, identifies clusters with metagenome-derived sequences.
#2) Outputs clusters containing ASVs and ASVs.
###############################################


#######
#Packages to be imported.
#######

import numpy as np
import Bio
from Bio import SeqIO
import re
import matplotlib.pyplot as plt
import csv
import pandas as pd

#######
#Data to be used.
#######

#MMseqs2 clustering results.
MMseq_data_path="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Clustering_by_mmseq\\V3_V4_ezbiocloud_plus_metagenome_Kaiju_derived_plus_Grotta_dataset_V3_V4_ASV_clust_res.tsv"
#ASV multifasta.
ASV_fasta_path="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Grotta_sponges_seqtab_nochim_mergers_primer_trimmed.fasta"
#ASV frequences.
ASV_freq_table_path="C:\Sutor\Science\Antarctic\Metagenome_based_16S_filtering\Grotta_sponges_seqtab_nochim_mergers.csv"

#Output:
Output_path="C:\Sutor\Science\Antarctic\\Metagenome_based_16S_filtering\\16S_filtering\\"


#######
#Read clustering results, prepare dictionary of clusters.
#######

def read_clusters_data(datapath):
    filein=open(datapath, 'r')
    Clusters_dict={}
    for line in filein:
        line=line.rstrip('\r\n').split('\t')
        cluster_rep=line[0]
        cluster_dep=line[1]
        #print(cluster_rep, cluster_dep)
        if cluster_rep in Clusters_dict:
            current_cluster_ar=Clusters_dict[cluster_rep]
            current_cluster_ar.append(cluster_dep)
            Clusters_dict[cluster_rep]=current_cluster_ar
            #print('Update cluster ', Clusters_dict[cluster_rep])
        else:
            Clusters_dict[cluster_rep]=[cluster_dep]
            #print('New cluster ', Clusters_dict[cluster_rep])
    filein.close()
    
    print('Number of clusters: ' + str(len(Clusters_dict)))
    return Clusters_dict


#######
#Analyse dictionary of clusters: clusters size distribution.
#######

def analyse_clusters_dict(Clusters_dict, output_path):
    Clusters_len=[]
    for cluster_rep, cluster_dep_ar in Clusters_dict.items():
        Clusters_len.append(len(cluster_dep_ar))
    
    
    #Plot distribution of clusters size.
    plt.hist(Clusters_len)
    plt.yscale('log')
    plt.xlabel('Size of clusters, 99% id')
    plt.ylabel('Number of clusters')
    plt.show()
    plt.savefig(output_path + 'MMseqs2_clusters_size_distribution_99id_Ezbiocloud_dp_plus_metagenome_derived_plus_Grotta_dataset.png')    

    return


#######
#Detect clusters containing metagenome-derived data: return clusters containing this data and clusters not containing separately.
#######

def detect_metagenome_data(Clusters_dict, pattern_to_detect):
    Seq_ids_matching_ar=[]
    Clusters_with_pattern={}
    Clusters_without_pattern={}
    for cluster_rep, cluster_dep_ar in Clusters_dict.items():
        pattern_search=0
        for seq_id in cluster_dep_ar:
            if re.match(pattern_to_detect, seq_id):
                Seq_ids_matching_ar.append(seq_id)
                pattern_search=1
        if pattern_search==0:
            Clusters_without_pattern[cluster_rep]=cluster_dep_ar
        else:
            Clusters_with_pattern[cluster_rep]=cluster_dep_ar
         
    print('Number of seq ids matching the pattern ' + str(pattern_to_detect) + ' : ' + str(len(Seq_ids_matching_ar)))           
    print('Number of clusters containing seq ids matching the pattern ' + str(pattern_to_detect) + ' : ' + str(len(Clusters_with_pattern)))           
    print('Number of clusters are not containing seq ids matching the pattern ' + str(pattern_to_detect) + ' : ' + str(len(Clusters_without_pattern)))
    return Clusters_with_pattern, Clusters_without_pattern


#######
#Return seq ids matching the pattern from clusters dictionary.
#Return representative sequence for clusters containing seq ids matching the pattern.
#######

def return_seq_ids_matching_pattern(Clusters_dict, pattern_to_detect):
    Seq_ids_matching_ar=[]
    Rep_seq_ar=[]
    for cluster_rep, cluster_dep_ar in Clusters_dict.items():
        pattern_search=0
        for seq_id in cluster_dep_ar:
            if re.match(pattern_to_detect, seq_id):
                Seq_ids_matching_ar.append(seq_id)
                pattern_search=1
        if pattern_search==1:
            Rep_seq_ar.append(cluster_rep)
            
    print('Number of seq ids matching the pattern ' + str(pattern_to_detect) + ' : ' + str(len(Seq_ids_matching_ar)))           
    print('Number of clusters containing seq ids matching the pattern ' + str(pattern_to_detect) + ' : ' + str(len(Rep_seq_ar)))
    return Seq_ids_matching_ar, Rep_seq_ar


#######
#Read initial ASV fasta.
#######

def read_mfa(ASV_fasta_path):
    filein=open(ASV_fasta_path, 'r')
    #Keep data in dictionary.
    ASV_dict={}
    for record in SeqIO.parse(filein, "fasta"):
        seq_16S=str(record.seq)
        id_16S=record.name
        ASV_dict[id_16S]=seq_16S
    filein.close()  
    return ASV_dict


#######
#Return sequence by ASV id, write final fasta file.
#######

def write_final_ASV_seqs(ASV_dict, Seq_ids_matching_ar, Output_path):
    fileout=open(Output_path+"Grotta_sponges_V3_V4_mergers_filtered_by_metagenome_data.fasta", 'w')
    i=0
    Seq_ids_matching_ar_ordered=[]
    for Seq_id, Seq_id_seq in ASV_dict.items():
        if Seq_id in Seq_ids_matching_ar:
            fileout.write(Seq_id + '\n' + Seq_id_seq + '\n')
            i+=1
            Seq_ids_matching_ar_ordered.append(Seq_id)
        else:
            print(Seq_id + ' was filtered out.')  
    fileout.close()
    print('Fraction of sequences survived after filtration: ' + str(float(len(Seq_ids_matching_ar))/len(ASV_dict)))
    print('Fraction of sequences survived after filtration and were written: ' + str(float(i)/len(ASV_dict)))
    return Seq_ids_matching_ar_ordered


#######
#Update table with ASV frequences, leave only survived ASV data.
#######

def update_ASV_freq_table(ASV_freq_table_path, Seq_ids_matching_ar_ordered, Output_path):
    ASV_freq_dataframe=pd.read_csv(ASV_freq_table_path, sep=',', header=0, index_col=0)
    #print(ASV_freq_dataframe)
    ASV_freq_dataframe_updated=ASV_freq_dataframe.loc[Seq_ids_matching_ar_ordered,:]
    #print(ASV_freq_dataframe_updated)   
    ASV_freq_dataframe_updated.to_csv(path_or_buf=Output_path+"Grotta_sponges_V3_V4_mergers_filtered_by_metagenome_data_frequences.csv", sep=',', quoting=csv.QUOTE_NONNUMERIC, header=True, index=True)
    return

MMseq_Clusters_dict=read_clusters_data(MMseq_data_path)
analyse_clusters_dict(MMseq_Clusters_dict, Output_path)
MMSeq_Clusters_with_pattern, MMSeq_Clusters_without_pattern=detect_metagenome_data(MMseq_Clusters_dict, "^\d+_\D+_\D+_\d+") #Searches for metagenome-derived seqs.
Seq_ids_matching_ar, Rep_seq_ar=return_seq_ids_matching_pattern(MMSeq_Clusters_without_pattern, '^ASV\d+') #Searches for ASVs.
ASV_dict=read_mfa(ASV_fasta_path)
Seq_ids_matching_ar_ordered=write_final_ASV_seqs(ASV_dict, Seq_ids_matching_ar, Output_path)
update_ASV_freq_table(ASV_freq_table_path, Seq_ids_matching_ar_ordered, Output_path)
