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
from collections import Counter
import csv
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.cluster.hierarchy as sch

################################################################################ 

# consolidateSeeds
# Consolidation seed weights into a single data frame.
# This function reads in the seed metabolites and their weights for each 
# genome. The data is read into a dataframe and then written to file. Structure:
# Rows: metabolites
# Columns: graphs
# Entries: unweighted seed set values

def consolidateSeeds(dirList, externalDataDir, processedDataDir, summaryStatsDir):
    
    print 'Consolidate seed sets'
        
# Create a data frame to store the data
    revEcolMatrixDF = pd.DataFrame(columns=["Metabolite"])

# Read in list of seed weights into a temporary data frame. Perform an outer
# join with the existing data frame to incorporate the new list of weights.
    for curDir in dirList:
        tempDF = pd.read_csv('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedWeights.txt', names=['Metabolite',curDir])
        revEcolMatrixDF = pd.merge(revEcolMatrixDF, tempDF, how='outer', on="Metabolite")

# Replace all the NaN values with zeros
    revEcolMatrixDF.fillna(0, inplace=True)

# Append a new column containing common names associated with metabolite IDs.
# The file metabMap.csv was created manually from the seed database, and should
# be updated to reflect the particulars of your data set.
    namesDF = pd.read_csv('../'+externalDataDir+'/'+'metabMap.csv', names=['Metabolite','CommonName'])
    revEcolMatrixDF = pd.merge(revEcolMatrixDF, namesDF, how='inner', on="Metabolite")

# Rearrange the order of the columns so that the common name is in front
    newOrder = revEcolMatrixDF.columns.tolist()
    newOrder = newOrder[-1:] + newOrder[:-1]
    revEcolMatrixDF = revEcolMatrixDF[newOrder]
    
# Export the matrix of seed weights
    pd.DataFrame.to_csv(revEcolMatrixDF, '../'+summaryStatsDir+'/'+'seedMatrixWeighted.csv')
    
    return revEcolMatrixDF

################################################################################ 

# normalizedSeedCounts
# This function reads in the seed compounds for all metabolic networks in the
# given directory, and computes the frequency with which a compound appears
# in the genomes.

def normalizedSeedCounts(dirList, processedDataDir, summaryStatsDir):

    numSubDir = len(dirList)
    print 'Computing normalized seed counts'

# Read in the seed sets from file
# seedSetList is a list of lists. Each outer list contains all the seed sets
# for that graph.
    seedSetList = []
    
    for curDir in dirList:

# read in list of seed sets for the current genome
        with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedSets.txt', 'r') as mySeedFile:
            reader = csv.reader(mySeedFile, delimiter=',')
            mySeedSet = list(reader)

# Append list of seeds to the master list
        seedSetList.append(mySeedSet)
    
# Count the number of occurences of each metabolite across all seed sets    
    seedMetabCount = Counter()
    for seedSet in seedSetList:
        for seed in seedSet:
            for metab in seed:
                seedMetabCount[metab] += 1
            
# Write a normalized list of counts to file. First normalize the counts. Then, 
# convert the Counter to a list of tuples (metabolite, counts) and write to file.
# Normalize
    for key, value in seedMetabCount.items():
        seedMetabCount[key] = float(value) / numSubDir

# Convert to a list of tuples
    seedMetabCountList = seedMetabCount.most_common()

# Write to file
    seedCounts = open('../'+summaryStatsDir+'/SeedCounts.txt', 'wb')
    seedCounts.write('Seed Metabolite,Frequency\n')
    writer = csv.writer(seedCounts, lineterminator='\n')
    writer.writerows(seedMetabCountList)
    seedCounts.close()
    
    return

################################################################################ 

# clusterSeedSets
# This function creates a dendrogram of the seed sets. The genomes are
# clustered based on weighted vectors of seed sets, using the euclidean
# distance and UPGMA (average linkage) clustering. A second dendrogram is
# constructed based on metabolite weights across genomes. Then, the matrix
# of seed weights is reordered to reflect the order of the dendrograms and is
# visualized. Genome names are colored as specified in colorFile.

def clusterSeedSets(revEcolMatrixDF, dirList, externalDataDir, summaryStatsDir, colorFile):

    print 'Computing dendrogram'
    
# Python clustering algorithms require the data to be an ndarray, with each
# row corresponding to a set of observations. Thus, calculation of the linkage
# for the genomes must be performed on the transpose.
    revEcolMatrix=pd.DataFrame.as_matrix(revEcolMatrixDF.drop(['CommonName','Metabolite'], 1))
    revEcolMatrixT=np.transpose(pd.DataFrame.as_matrix(revEcolMatrixDF.drop(['CommonName','Metabolite'], 1)))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7.5, 61))

# Compute and plot dendrogram for genomes, which will be above the graph.
# Define the size of the dendrogram
    ax1 = fig.add_axes([0.05,0.9,0.8,0.05], frame_on=False)
# Compute the linkage matrix
    genomeLinkage = sch.linkage(revEcolMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    genomeClust = sch.dendrogram(genomeLinkage, distance_sort='True', color_threshold=0)
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Compute and plot dendrogram for metabolites which will be beside the graph. 
# Define the size of the dendrogram
    ax2 = fig.add_axes([0,0.1,0.05,0.8], frame_on=False)
# Compute the linkage matrix
    metabLinkage = sch.linkage(revEcolMatrix, method='average', metric='euclidean')
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
    revEcolMatrix = revEcolMatrix[:,idx1]
    revEcolMatrix = revEcolMatrix[idx2,:]
# Plot the weight matrix
    im = axmatrix.matshow(revEcolMatrix, aspect='auto', origin='lower')
# No tick marks along axes
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+colorFile)
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]

# Add genome names to the bottom axis
    axmatrix.set_xticks(range(len(revEcolMatrixT)))
    axmatrix.set_xticklabels([dirList[i] for i in idx1], minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(axmatrix.get_xticklabels(), genomeColors):
        xtick.set_color(color)

# Add metabolite names to the right axis
    axmatrix.set_yticks(range(len(revEcolMatrix)))
    axmatrix.set_yticklabels([revEcolMatrixDF['CommonName'][i] for i in idx2], minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    plt.yticks(fontsize=8)

# Plot colorbar.
    axcolor = fig.add_axes([0.97,0.9,0.03,0.05])
    plt.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('../'+summaryStatsDir+'/'+'seedSetDendrogram.png')

################################################################################ 

# clusterOnly
# The dendrogram of the seed sets is quite large, making it difficult to see
# how the genome cluster. This function displays just the dendrogram of the 
# genomes, based on weighted vectors of seed sets.

def clusterOnly(revEcolMatrixDF, dirList, externalDataDir, colorFile):

    revEcolMatrix=pd.DataFrame.as_matrix(revEcolMatrixDF.drop(['CommonName','Metabolite'], 1))
    revEcolMatrixT=np.transpose(pd.DataFrame.as_matrix(revEcolMatrixDF.drop(['CommonName','Metabolite'], 1)))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7,2))
    
# Define the size of the dendrogram
    ax1 = fig.add_axes([0, 0, 1, 1], frame_on=False)
# Compute the linkage matrix
    genomeLinkage = sch.linkage(revEcolMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    genomeClust = sch.dendrogram(genomeLinkage, distance_sort='True', color_threshold=0)
    idx1 = genomeClust['leaves']
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+colorFile)
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]
    
    axisLength = 1- (1/float(len(dirList)))
    axisStart = (1 - axisLength) / 2

# Add genome names to the bottom axis
    ax2 = fig.add_axes([axisStart, 0, axisLength, 0], frame_on=False)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xticks(range(len(revEcolMatrixT)))
    ax2.set_xticklabels([dirList[i] for i in idx1], minor=False)
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(ax2.get_xticklabels(), genomeColors):
        xtick.set_color(color)

################################################################################ 

# Reverse Ecology Metric: Metabolic Competition Index
# computeMetabCompete

# The metabolic competition index (MII) indicates the competitive potential 
# between two organisms. For two organisms A and B, MII is calculated as the 
# fraction of compounds in the seed set of A which are also in seed set of B,
# and ranges from 0 to 1. Because seed compounds are associated with a weight,
# this metric is calculated as a normalized weighted sum. Note that this 
# metric is non-symmetric.

# This function computes the MII between all pairs of genomes. A nested loop
# is used to iterate all over all pairs of genomes, and a dataframe is used to
# store the MII for each pair. First, the weighted seed sets for each genome
# are read into dataframes. The overlap between sets is performed using an
# inner join on the two seed sets. The inner join contains seed weights for
# both genomes; weights from the first genome are used when computing the MII.
# The MII is then computed as the ratio of two weighted sums: the seed set
# overlap between A and B, and the seed set of A. After all pairwise scores 
# are computed, the output is written to file.

def computeMetabCompete(dirList, processedDataDir, summaryStatsDir):
    
    numSubDir = len(dirList)
# Establish a matrix to store MII scores
    metabCompete = pd.DataFrame(np.zeros((numSubDir, numSubDir)), index=dirList, columns=dirList)

# Use a nested loop to loop over all genome pairs. outerDir coresponds to
# organism A and innerDir to organism B
    for outerDir in dirList:
        for innerDir in dirList:

# Read in the list of seed sets and their weights for organisms A and B
            seedWeightOuter = pd.read_csv('../'+processedDataDir+'/'+outerDir+'/'+outerDir+'SeedWeights.txt', 
                                      header=None, names=['Metabolite', 'Outer Weight'])
            seedWeightInner = pd.read_csv('../'+processedDataDir+'/'+innerDir+'/'+innerDir+'SeedWeights.txt', 
                                      header=None, names=['Metabolite', 'Inner Weight'])

# Compute the overlap between seed sets using an inner join
            overlapSeeds = pd.merge(seedWeightOuter, seedWeightInner, on='Metabolite')

# Sum seed compound weights for the overlap between A and B
            upperSum = overlapSeeds.loc[:,'Outer Weight'].sum()

# Sum seed compound weights for the seed set of A
            lowerSum = seedWeightOuter.loc[:,'Outer Weight'].sum()

# The MII is the ratio of these two values. Compute and store.
            metabCompete.loc[outerDir, innerDir] = upperSum / lowerSum

# When loop complete, write to file
    metabCompete.to_csv('../'+summaryStatsDir+'/'+'MetabolicCompetitionScores.csv')  
    
    return metabCompete

################################################################################ 

# Reverse Ecology Metric: Metabolic Complementarity Index
# computeMetabComplement

# The metabolic complementarity index (MCI) indicates the complementarity of
# two organisms' niches. For two organisms A and B, MCI is calculated as the 
# fraction of compounds in the seed set of A which are in the metabolic network
# of B, but not in B's seed set. MCI can range from 0 to 1. Note that this 
# metric is non-symmetric.

# This function computes the MCI between all pairs of genomes. A nested loop
# is used to iterate all over all pairs of genomes, and a dataframe is used to
# store the MCI for each pair. First, the weighted seed sets for each genome
# are read into dataframes. Then, the metabolic network of B is read in as a 
# digraph and convert to a dataframe. The list of non-seed compounds for
# organism B is computed from these two dataframes. Then, the overlap between 
# A's seed compounds and B's non-seed compounds is performed using an inner 
# join on the two seed sets. The MCI is then computed as the ratio of two
# sizes: the size of A's seed set, and the size of the overlap with B's non-
# seeds. After all pairwise scores are computed, the output is written to file.

def computeMetabComplement(dirList, processedDataDir, summaryStatsDir):
    
    numSubDir = len(dirList)
    
# Establish a matrix to store MCI scores
    metabComplement = pd.DataFrame(np.zeros((numSubDir, numSubDir)), index=dirList, columns=dirList)

# Use a nested loop to loop over all genome pairs. outerDir coresponds to
# organism A and innerDir to organism B
    for outerDir in dirList:
        for innerDir in dirList:

# Read in the list of seed sets and their weights for organisms A and B
            seedWeightOuter = pd.read_csv('../'+processedDataDir+'/'+outerDir+'/'+outerDir+'SeedWeights.txt', header=None, names=['Metabolite', 'Outer Weight'])

# Read in the lists of seed compounds and all compounds for the inner genome
            seedWeightInner = pd.read_csv('../'+processedDataDir+'/'+innerDir+'/'+innerDir+'SeedWeights.txt', header=None, names=['Metabolite', 'Inner Weight'])

# Read in the metabolic network of B via its adjancency list representation and
# converts it to a pandas Dataframe. First, the list of nodes is extracted from
# the graph. The list is then converted to a pandas Series, embedded in a dict,
# and used to create a dataframe of the nodes.
            allNodesInner = pd.DataFrame({'Metabolite' : pd.Series(nx.read_adjlist('../'+processedDataDir+'/'+innerDir+'/'+innerDir+'AdjList.txt',delimiter='\t', create_using=nx.DiGraph()).nodes())})

# Compute the list of non-seed compounds for organisms B
            nonSeedsInner = allNodesInner[~allNodesInner.Metabolite.isin(seedWeightInner.Metabolite)]

# Compute the overlap between A's seeds and B's non-seeds
            overlapSeeds = pd.merge(seedWeightOuter, nonSeedsInner, on='Metabolite')

# Compute the ratio of these two sets
            metabComplement.loc[outerDir, innerDir] = float(len(overlapSeeds)) / float(len(seedWeightOuter))

# When loop complete, write to file
    metabComplement.to_csv('../'+summaryStatsDir+'/'+'MetabolicComplementarityScores.csv')
    
    return metabComplement
    
################################################################################ 
    
# Clustering and Visualization of a matrix dataframe
# clusterPairwise

# This function creates a dendrogram of pairwise reverse ecology metric, such 
# as metabolic complementarity scores. The genomes are clustered using the 
# euclidean distance and UPGMA (average linkage) clustering. Because the scores 
# are non-symmetric, each axis will have separate clustering. Genome names are
# colored as specified in colorFile.

def clusterPairwise(revEcolMatrixDF, dirList, externalDataDir, summaryStatsDir, colorFile, fileName):
# Python clustering algorithms require the data to be an ndarray, with each
# row corresponding to a set of observations.
    revEcolMatrix=pd.DataFrame.as_matrix(revEcolMatrixDF)
    revEcolMatrixT=np.transpose(pd.DataFrame.as_matrix(revEcolMatrixDF))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7.5, 7.5))

# Compute and plot the first dendrogram, which will be above the graph.
# Define the size of the dendrogram
    ax1 = fig.add_axes([0.1, 0.9, 0.8, 0.1], frame_on=False)
# Compute the linkage matrix
    firstLinkage = sch.linkage(revEcolMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    firstClust = sch.dendrogram(firstLinkage, distance_sort='True', color_threshold=0)
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Compute and plot the second dendrogram, which will be beside the graph.
# Define the size of the dendrogram
    ax2 = fig.add_axes([0, 0.1, 0.1, 0.8], frame_on=False)
# Compute the linkage matrix
    secondLinkage = sch.linkage(revEcolMatrix, method='average', metric='euclidean')
# Compute the dendrogram
    secondClust = sch.dendrogram(secondLinkage, orientation='right', color_threshold=0)
# No tick marks along axes
    ax2.set_xticks([])
    ax2.set_yticks([])

# Plot the matrix of competition scores
# Define the size of the plot
    axmatrix = fig.add_axes([0.1, 0.1, 0.8, 0.8])
# Determine ordering of the dendrograms and rearrange the weight matrix
    idx1 = firstClust['leaves']
    idx2 = secondClust['leaves']
    revEcolMatrix = revEcolMatrix[:,idx1]
    revEcolMatrix = revEcolMatrix[idx2,:]
# Plot the weight matrix
    im = axmatrix.matshow(revEcolMatrix, aspect='auto', origin='lower')
# No tick marks along axes
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+colorFile)
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]

# Add genome names to the bottom axis
    axmatrix.set_xticks(range(len(revEcolMatrixT)))
    axmatrix.set_xticklabels([dirList[i] for i in idx1], minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(axmatrix.get_xticklabels(), genomeColors):
        xtick.set_color(color)

# Add genome names to the right axis
    axmatrix.set_yticks(range(len(revEcolMatrix)))
    axmatrix.set_yticklabels([dirList[i] for i in idx2], minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    plt.yticks(fontsize=8)
    for ytick, color in zip(axmatrix.get_yticklabels(), genomeColors):
        ytick.set_color(color)

# Plot colorbar.
    axcolor = fig.add_axes([0.9, 0.9, 0.03, 0.1])
    plt.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('../'+summaryStatsDir+'/'+fileName)

    return