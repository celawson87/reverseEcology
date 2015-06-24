###############################################################################
# seedFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with seed sets.
################################################################################ 

# Import Python packages.
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch


#%% Consolidation seed weights into a single data frame.
# This code snippet reads in the seed metabolites and their weights for each 
# genome. The data is read into a dataframe and then written to file. Structure:
# Rows: metabolites
# Columns: graphs
# Entries: unweighted seed set values

def consolidateSeeds(dirList, externalDataDir, processedDataDir, summaryStatsDir):
    
    print 'Consolidate seed sets'
        
# Create a data frame to store the data
    seedMatrixDF = pd.DataFrame(columns=["Metabolite"])

# Read in list of seed weights into a temporary data frame. Perform an outer
# join with the existing data frame to incorporate the new list of weights.
    for curDir in dirList:
        tempDF = pd.read_csv('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedWeights.txt', names=['Metabolite',curDir])
        seedMatrixDF = pd.merge(seedMatrixDF, tempDF, how='outer', on="Metabolite")

# Replace all the NaN values with zeros
    seedMatrixDF.fillna(0, inplace=True)

# Append a new column containing common names associated with metabolite IDs.
# The file metabMap.csv was created manually from the seed database, and should
# be updated to reflect the particulars of your data set.
    namesDF = pd.read_csv('../'+externalDataDir+'/'+'metabMap.csv', names=['Metabolite','CommonName'])
    seedMatrixDF = pd.merge(seedMatrixDF, namesDF, how='inner', on="Metabolite")

# Rearrange the order of the columns so that the common name is in front
    newOrder = seedMatrixDF.columns.tolist()
    newOrder = newOrder[-1:] + newOrder[:-1]
    seedMatrixDF = seedMatrixDF[newOrder]
    
# Export the matrix of seed weights
    pd.DataFrame.to_csv(seedMatrixDF, '../'+summaryStatsDir+'/'+'seedMatrixWeighted.csv')
    
    return seedMatrixDF


#%% Clustering and Visualization of Seed Sets

def clusterSeedSets(seedMatrixDF, dirList, externalDataDir, summaryStatsDir):
# This code snippet creates a dendrogram of the seed sets. The genomes are
# clustered based on weighted vectors of seed sets, using the euclidean
# distance and UPGMA (average linkage) clustering. A second dendrogram is
# constructed based on metabolite weights across genomes. Then, the matrix
# of seed weights is reordered to reflect the order of the dendrograms and is
# visualized. Genome names are colored according to their lineage, with acI
# sub-divided into acI-A and acI-B.

    print 'Computing dendrogram'
    
# Python clustering algorithms require the data to be an ndarray, with each
# row corresponding to a set of observations. Thus, calculation of the linkage
# for the genomes must be performed on the transpose.
    seedMatrix=pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1))
    seedMatrixT=np.transpose(pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1)))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7.5, 61))

# Compute and plot dendrogram for genomes, which will be above the graph.
# Define the size of the dendrogram
    ax1 = fig.add_axes([0.05,0.9,0.8,0.05], frame_on=False)
# Compute the linkage matrix
    genomeLinkage = sch.linkage(seedMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    genomeClust = sch.dendrogram(genomeLinkage, distance_sort='True', color_threshold=0)
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Compute and plot dendrogram for metabolites which will be beside the graph. 
# Define the size of the dendrogram
    ax2 = fig.add_axes([0,0.1,0.05,0.8], frame_on=False)
# Compute the linkage matrix
    metabLinkage = sch.linkage(seedMatrix, method='average', metric='euclidean')
# Compute the dendrogram
    metabClust = sch.dendrogram(metabLinkage, orientation='right', color_threshold=0)
# No tick marks along axes
    ax2.set_xticks([])
    ax2.set_yticks([])
    
# Plot the matrix of seed weights.
# Define the size of the plot
    axmatrix = fig.add_axes([0.05,0.1,0.8,0.8])
# Determine ordering of the dendrograms and rearrange the weight matrix
    idx1 = genomeClust['leaves']
    idx2 = metabClust['leaves']
    seedMatrix = seedMatrix[:,idx1]
    seedMatrix = seedMatrix[idx2,:]
# Plot the weight matrix
    im = axmatrix.matshow(seedMatrix, aspect='auto', origin='lower')
# No tick marks along axes
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+'tribalColors.csv')
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]

# Add genome names to the bottom axis
    axmatrix.set_xticks(range(len(seedMatrixT)))
    axmatrix.set_xticklabels([dirList[i] for i in idx1], minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(axmatrix.get_xticklabels(), genomeColors):
        xtick.set_color(color)

# Add metabolite names to the right axis
    axmatrix.set_yticks(range(len(seedMatrix)))
    axmatrix.set_yticklabels([seedMatrixDF['CommonName'][i] for i in idx2], minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    plt.yticks(fontsize=8)

# Plot colorbar.
    axcolor = fig.add_axes([0.97,0.9,0.03,0.05])
    plt.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('../'+summaryStatsDir+'/'+'seedSetDendrogram.png')

#%% Clustering and Visualization of Seed Sets

def clusterOnly(seedMatrixDF, dirList, externalDataDir):

    seedMatrix=pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1))
    seedMatrixT=np.transpose(pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1)))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7,2))
    
# Define the size of the dendrogram
    ax1 = fig.add_axes([0, 0, 1, 1], frame_on=False)
# Compute the linkage matrix
    genomeLinkage = sch.linkage(seedMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    genomeClust = sch.dendrogram(genomeLinkage, distance_sort='True', color_threshold=0)
    idx1 = genomeClust['leaves']
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+'tribalColors.csv')
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]
    
    axisLength = 1- (1/float(len(dirList)))
    axisStart = (1 - axisLength) / 2
# Add genome names to the bottom axis
    ax2 = fig.add_axes([axisStart, 0, axisLength, 0], frame_on=False)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xticks(range(len(seedMatrixT)))
    ax2.set_xticklabels([dirList[i] for i in idx1], minor=False)
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(ax2.get_xticklabels(), genomeColors):
        xtick.set_color(color)
