###############################################
##Dmitry Sutormin, 2020##
##Antarctic metagenome analysis##

#Takes table with taxons abundance for different steps of a pipeline and visualize the abundance.
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


#Path to taxa table.
Taxa_table="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Antarctic\Top_genus\Top_genus_for_different_stages_extended.xlsx"

No_pooling=pd.read_excel(Taxa_table, sheet_name='No_pooling', header=0)
With_pooling=pd.read_excel(Taxa_table, sheet_name='With_pooling', header=0)

###############
#Create list of random colors.
###############

def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    ###Stolen from here: https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np
    
    if type not in ('bright', 'soft'):
        print('Please choose "bright" or "soft" for type')
        return
    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type=='bright':
        randRGBcolors=[(np.random.uniform(low=0.0, high=1), np.random.uniform(low=0, high=1), np.random.uniform(low=0, high=1)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0]=[0, 0, 0]

        if last_color_black:
            randRGBcolors[-1]=[0, 0, 0]

        random_colormap=LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type=='soft':
        low=0.6
        high=0.95
        randRGBcolors=[(np.random.uniform(low=low, high=high), np.random.uniform(low=low, high=high), np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0]=[0, 0, 0]
        
        if last_color_black:
            randRGBcolors[-1]=[0, 0, 0]
        random_colormap=LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)
    
        # Display colorbar
        if verbose:
            from matplotlib import colors, colorbar
            from matplotlib import pyplot as plt
            fig, ax=plt.subplots(1, 1, figsize=(15, 0.5))
            bounds=np.linspace(0, nlabels, nlabels + 1)
            norm=colors.BoundaryNorm(bounds, nlabels)
            cb=colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                       boundaries=bounds, format='%1i', orientation=u'horizontal')
    return randRGBcolors    

##############
#Prepare list of items (taxons) and assign them colors.
##############

def colors_assignment(names_list_1, names_list_2):
    #Create dictionary containing all names from both lists.
    names_list_combined={}
    for name in names_list_1:
        names_list_combined[name]=0
    for name in names_list_2:
        if name not in names_list_combined:
            names_list_combined[name]=0
    #Create list of colors.
    colors=rand_cmap(len(names_list_combined), type='bright', first_color_black=False, last_color_black=True, verbose=False)
    #Assign colors to dictionary keys.
    i=0
    for key in names_list_combined:
        names_list_combined[key]=colors[i]
        i+=1
    return names_list_combined

#Plot data.
fig, plot_av=plt.subplots(1,1,figsize=(12,6), dpi=100)
#Prepare colors.
names_and_colors=colors_assignment(No_pooling['Genus'].tolist(), With_pooling['Genus'].tolist())
#Prepare points size.
No_pooling_Before_MMSeqs_size=np.sqrt(np.array(No_pooling['Before MMSeqs'].tolist())/np.pi)
No_pooling_After_MMSeqs_size=np.sqrt(np.array(No_pooling['After MMSeqs'].tolist())/np.pi)
No_pooling_After_Decontam_size=np.sqrt(np.array(No_pooling['After Decontam'].tolist())/np.pi)


#Place points and legend no pooling.
points_1_ar=[]
for i in range(len(No_pooling_Before_MMSeqs_size)):
    #Points.
    plot_av.scatter(1, No_pooling['Before MMSeqs'].tolist()[i], s=No_pooling_Before_MMSeqs_size[i], color=names_and_colors[No_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.2, alpha=1, zorder=1000)
    plot_av.scatter(1.6, No_pooling['After MMSeqs'].tolist()[i], s=No_pooling_After_MMSeqs_size[i], color=names_and_colors[No_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.2, alpha=1, zorder=1000)
    plot_av.scatter(2.2, No_pooling['After Decontam'].tolist()[i], s=No_pooling_After_Decontam_size[i], color=names_and_colors[No_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.2, alpha=1, zorder=1000)
    #Points annotation.
    if i%2==0:
        print(i, 'even', No_pooling['Genus'].tolist()[i])
        plot_av.annotate(No_pooling['Genus'].tolist()[i], xy=(2.25, No_pooling['After Decontam'].tolist()[i]), xycoords='data', va="center", fontsize=2)
    elif i%2==1:
        print(i, 'odd', No_pooling['Genus'].tolist()[i])
        plot_av.annotate(No_pooling['Genus'].tolist()[i], xy=(2.38, No_pooling['After Decontam'].tolist()[i]), xycoords='data', va="center", fontsize=2)
    else:
        print(i, 'else', No_pooling['Genus'].tolist()[i])   
    #Legend.
    points_1=plot_av.scatter(-1, No_pooling['After Decontam'].tolist()[i], s=80, color=names_and_colors[No_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.4, zorder=1000, label=No_pooling['Genus'].tolist()[i])
    points_1_ar.append(points_1)
#Place links no pooling.
for i in range(len(No_pooling_Before_MMSeqs_size)):
    if No_pooling_Before_MMSeqs_size[i]!=0 and No_pooling_After_MMSeqs_size[i]!=0:
        plot_av.plot([1,1.6], [No_pooling['Before MMSeqs'].tolist()[i], No_pooling['After MMSeqs'].tolist()[i]], 'k-', linewidth=0.3, zorder=10)
    if No_pooling_After_MMSeqs_size[i]!=0 and No_pooling_After_Decontam_size[i]!=0:
        plot_av.plot([1.6,2.2], [No_pooling['After MMSeqs'].tolist()[i], No_pooling['After Decontam'].tolist()[i]], 'k-', linewidth=0.3, zorder=10)    


#Plot data pooling.
#Prepare points size.
Pooling_Before_MMSeqs_size=np.sqrt(np.array(With_pooling['Before MMSeqs'].tolist())/np.pi)
Pooling_After_MMSeqs_size=np.sqrt(np.array(With_pooling['After MMSeqs'].tolist())/np.pi)
Pooling_After_Decontam_size=np.sqrt(np.array(With_pooling['After Decontam'].tolist())/np.pi)

#Place points and legend pooling.
points_2_ar=[]
for i in range(len(Pooling_Before_MMSeqs_size)):
    #Points.
    plot_av.scatter(2.8, With_pooling['Before MMSeqs'].tolist()[i], s=Pooling_Before_MMSeqs_size[i], color=names_and_colors[With_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.2, alpha=1, zorder=1000)
    plot_av.scatter(3.4, With_pooling['After MMSeqs'].tolist()[i], s=Pooling_After_MMSeqs_size[i], color=names_and_colors[With_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.2, alpha=1, zorder=1000)
    plot_av.scatter(4.0, With_pooling['After Decontam'].tolist()[i], s=Pooling_After_Decontam_size[i], color=names_and_colors[With_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.2, alpha=1, zorder=1000)
    #Points annotation.
    if i%2==0:
        print(i, 'even', With_pooling['Genus'].tolist()[i])
        plot_av.annotate(With_pooling['Genus'].tolist()[i], xy=(4.05, With_pooling['After Decontam'].tolist()[i]), xycoords='data', va="center", fontsize=2)
    elif i%2==1:
        print(i, 'odd', With_pooling['Genus'].tolist()[i])
        plot_av.annotate(With_pooling['Genus'].tolist()[i], xy=(4.18, With_pooling['After Decontam'].tolist()[i]), xycoords='data', va="center", fontsize=2)
    else:
        print(i, 'else', With_pooling['Genus'].tolist()[i])    
    #Legend.
    points_2=plot_av.scatter(-1, With_pooling['After Decontam'].tolist()[i], s=80, color=names_and_colors[With_pooling['Genus'].tolist()[i]], edgecolors='black', linewidth=0.4, zorder=1000, label=With_pooling['Genus'].tolist()[i])
    points_2_ar.append(points_2)
#Place links no pooling.
for i in range(len(Pooling_Before_MMSeqs_size)):
    if Pooling_Before_MMSeqs_size[i]!=0 and Pooling_After_MMSeqs_size[i]!=0:
        plot_av.plot([2.8,3.4], [With_pooling['Before MMSeqs'].tolist()[i], With_pooling['After MMSeqs'].tolist()[i]], 'k-', linewidth=0.3, zorder=10)
    if Pooling_After_MMSeqs_size[i]!=0 and Pooling_After_Decontam_size[i]!=0:
        plot_av.plot([3.4,4.0], [With_pooling['After MMSeqs'].tolist()[i], With_pooling['After Decontam'].tolist()[i]], 'k-', linewidth=0.3, zorder=10)    

plot_av.set_ylabel('Abundance', size=20)
plot_av.set_xticks([1, 1.6, 2.2, 2.8, 3.4, 4.0], minor=False)
plot_av.set_xticklabels(['Before MMSeqs\nNo pooling', 'After MMSeqs\nNo pooling', 'After Decontam\nNo pooling', 'Before MMSeqs\nWith pooling', 'After MMSeqs\nWith pooling', 'After Decontam\nWith pooling'], rotation=0, size=12)
plot_av.set_ylim([1, 5*max(max(With_pooling['Before MMSeqs']), max(With_pooling['After MMSeqs']), max(With_pooling['After Decontam']),
                              max(No_pooling['Before MMSeqs']), max(No_pooling['After MMSeqs']), max(No_pooling['After Decontam']))])
plot_av.set_xlim([0.8,4.3])
plot_av.set_yscale('log')
#Configure legend.
lab1 = [ h.get_label() for h in points_1_ar]
lab2 = [ h.get_label() for h in points_2_ar]
leg11 = plt.legend(points_1_ar[:int(len(points_1_ar)/3)], lab1[:int(len(points_1_ar)/3)], loc='upper center', fontsize=6, bbox_to_anchor=(0.06,0.3))
leg12 = plt.legend(points_1_ar[int(len(points_1_ar)/3):int(2*len(points_1_ar)/3)], lab1[int(len(points_1_ar)/3):int(2*len(points_1_ar)/3)], loc='upper center', fontsize=6, bbox_to_anchor=(0.17,0.3))
leg13 = plt.legend(points_1_ar[int(2*len(points_1_ar)/3):], lab1[int(2*len(points_1_ar)/3):], loc='upper center', fontsize=6, bbox_to_anchor=(0.28,0.3))

#leg2 = plt.legend(points_2_ar, lab2, loc='upper right', fontsize=6, bbox_to_anchor=(0.99,0.99))
leg21 = plt.legend(points_2_ar[:int(len(points_2_ar)/4)+1], lab2[:int(len(points_2_ar)/4)+1], loc='upper center', fontsize=6, bbox_to_anchor=(0.07,0.325))
leg22 = plt.legend(points_2_ar[int(len(points_2_ar)/4)+1:int(2*len(points_2_ar)/4)+1], lab2[int(len(points_2_ar)/4)+1:int(2*len(points_2_ar)/4)+1], loc='upper center', fontsize=6, bbox_to_anchor=(0.25,0.325))
leg23 = plt.legend(points_2_ar[int(2*len(points_2_ar)/4)+1:int(3*len(points_2_ar)/4)+2], lab2[int(2*len(points_2_ar)/4)+1:int(3*len(points_2_ar)/4)+2], loc='upper center', fontsize=6, bbox_to_anchor=(0.43,0.325))
leg24 = plt.legend(points_2_ar[int(3*len(points_2_ar)/4)+2:], lab2[int(3*len(points_2_ar)/4)+2:], loc='upper center', fontsize=6, bbox_to_anchor=(0.55,0.325))

plt.gca().add_artist(leg21)
plt.gca().add_artist(leg22)
plt.gca().add_artist(leg23)
plt.gca().add_artist(leg24)
#plt.legend()
plt.tight_layout()
plt.show()
plt.savefig('C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\Antarctic\Top_genus\Taxa_filtration_abundance_extended.png', dpi=300, size=(12,13))

#colors=rand_cmap(len(No_pooling_Before_MMSeqs_size), type='bright', first_color_black=False, last_color_black=True, verbose=False)