###############################################
##Dmitry Sutormin, 2020##
##Antarctic metagenome analysis##

#Takes table with taxons abundance for different sampling points and performs clustering.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import binom
import scipy.cluster.hierarchy as sch

#Path to taxa table.
Taxa_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Antarctic\Top_genus\Taxonomy_tables\phyl_pool_after_Decontam.xlsx"

Abundance_df=pd.read_excel(Taxa_table, sheet_name='Abundance', header=0, index_col=0)


#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(clust_matrix, outpath_folder, file_name):
    X = clust_matrix.values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [clust_matrix.columns.tolist()[i] for i in list((np.argsort(ind)))]
    clust_matrix = clust_matrix.reindex(columns, axis=1)
    clust_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    return clust_matrix

#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def make_correlation_matrix_plot(df, cor_method, title, outpath_folder, file_name):
    fig=plt.figure(figsize=(10,10), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #Create correlation matrix and heatmap.
    df_cor_matrix=df.corr(method=cor_method)
    df_cor_matrix.to_csv(outpath_folder+file_name+'.csv', sep='\t', header=True, index=True)
    color_ax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1, aspect="equal")
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    #Label ticks.
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    ax1.set_ylim(sorted(ax1.get_xlim(), reverse=True)) #Solves a bug in matplotlib 3.1.1 discussed here: https://stackoverflow.com/questions/56942670/matplotlib-seaborn-first-and-last-row-cut-in-half-of-heatmap-plot
    #Create text annotation for heatmap pixels.
    #for i in range(len(labels)):
    #    for j in range(len(labels)):
    #        text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 2), ha="center", va="center", color="black")    
    #Add colorbar.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    #axins1=inset_axes(ax1, width="5%",  height="50%",  loc='upper right', bbox_to_anchor=(1.05, 0., 1, 1), bbox_transform=ax1.transAxes, borderpad=0)    #From here: https://matplotlib.org/3.1.1/gallery/axes_grid1/demo_colorbar_with_inset_locator.html 
    fig.colorbar(color_ax, ticks=[-1.00, -0.80, -0.60, -0.40, -0.20, 0.00, 0.20, 0.40, 0.60, 0.80, 1.00], shrink=0.5, panchor=(1.0, 0.5))
    plt.tight_layout()
    plt.savefig(outpath_folder+file_name+'.png', dpi=400, figsize=(10, 10))
    plt.show()
    plt.close()
    return df_cor_matrix

make_correlation_matrix_plot(Abundance_df, 'pearson', 'Sampling points', 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Antarctic\Top_genus\Taxonomy_tables\\', 'Sampling_points_after_Decontam_with_pooling_correlation_pearson')

Abundance_df_clusterized=Clustering(Abundance_df, 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Antarctic\Top_genus\Taxonomy_tables\\', 'Sampling_points_clusterized_after_Decontam_with_pooling_correlation_pearson')

make_correlation_matrix_plot(Abundance_df_clusterized, 'pearson', 'Sampling points clusterized', 'C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Antarctic\Top_genus\Taxonomy_tables\\', 'Sampling_points_clusterized_after_Decontam_with_pooling_correlation_pearson')