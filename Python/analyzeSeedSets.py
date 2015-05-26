###############################################################################
# analyzeSeedSets.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Script to import seed sets for each model and generate some basic plots.
# Plots include:
#  Count occurence of each metabolite across all seed sets
#  Generate clustergram of seed sets
################################################################################

#%% Preliminaries

# Import Python packages.
import csv
import os
from collections import Counter
from matplotlib import pyplot
import numpy as np

# Define data directories
# The user should ensure these directories match their data organization scheme
# processedDataDir: directory containing a folder for each genome. The folder
# name should match the genome name, and should contain a genome-scale model
# of the genome. The model should be named 'genomeNameBalanced.xml'.
# summaryStatsDir: directory for storing summary statistics
processedDataDir = 'ProcessedModelFiles'
summaryStatsDir = 'DataSummaries'

# Retrieve list of genomes to process by examing the contents of 
# 'processedDataDir', ignoring hidden folders. Subsequent computations are 
# performed by iterating over this list.
dirList =[]
for item in os.listdir('../'+processedDataDir):
    if not item.startswith('.'):
        dirList.append(item)

numSubDir = len(dirList)

#%% Read in the seed sets from file

# seedSetList is a list of lists. Each outer list contains all the seed sets
# for that graph.
seedSetList = []

# Read in the list of seed metabolites for each genome.
for curDir in dirList:
# read in list of seed sets for the current genome
    with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedSets.txt', 'r') as mySeedFile:
        reader = csv.reader(mySeedFile, delimiter=',')
        mySeedSet = list(reader)
# Append list of seeds to the master list
    seedSetList.append(mySeedSet)
    
#%% Calculate and plot a frequency distribubtion of metabolite occurances
# across seed sets.
    
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

# Turns out a plot of (name, frequency) pairs isn't very useful because
# there are so many of them!
#x = np.linspace(1, len(seedMetabCount), len(seedMetabCount))
#y = zip(*seedMetabCountList)[1]
#labels = zip(*seedMetabCountList)[0]

#pyplot.figure(figsize=(20,10))
#pyplot.scatter(x, y, marker='.')
#pyplot.xlim(0, len(x))
#pyplot.xticks(x, labels, rotation=90)
#pyplot.ylim(0, max(y))
#pyplot.xlabel('Metabolite')
#pyplot.ylabel('Frequency as a Seed')


#%% Clustering.

# A data frame would be a great way to visualize the data. Then we can extract
# the columns for clustering.
# Rows: metabolites
# Columns: graphs
# Entries: (un)weighted seed set values

# Import various packages
from pandas import DataFrame, read_csv
import pandas

# Create a data frame to store the data
seedMatrixDF = pandas.DataFrame(columns=["Metabolite"])

# Read in list of seed weights into a temporary data frame. Perform an outer
# join with the existing data frame to incorporate the new list of weights.
for curDir in dirList:
    tempDF = pandas.read_csv('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedWeights.txt', names=['Metabolite',curDir])
    seedMatrixDF = pandas.merge(seedMatrixDF, tempDF, how='outer', on="Metabolite")

# Replace all the NaN values with zeros
seedMatrixDF.fillna(0, inplace=True)


#%% Clustering and Visualization

from scipy import spatial
from scipy import cluster
import scipy
import scipy.cluster.hierarchy as sch

# Export the matrix of seed weights
pandas.DataFrame.to_csv(seedMatrixDF, '../'+summaryStatsDir+'/'+'seedMatrixWeighted.csv')

# Python clustering algorithms require the data to be an ndarray, with each
# row corresponding to a set of observations. Thus, calculation of the linkage
# for the genomes must be performed on the transpose.
seedMatrix=DataFrame.as_matrix(seedMatrixDF.drop('Metabolite', 1))
seedMatrixT=np.transpose(DataFrame.as_matrix(seedMatrixDF.drop('Metabolite', 1)))

#%% Clustering.

fig = pyplot.figure(figsize=(7.5, 61))

# Compute and plot dendogram for genomes, which will be above the graph. This
# dendogram is for the genomes.
ax1 = fig.add_axes([0.05,0.9,0.8,0.05], frame_on=False)
Y = sch.linkage(seedMatrixT, method='complete', metric='euclidean')
Z1 = sch.dendrogram(Y, distance_sort='True')
ax1.set_xticks([])
ax1.set_yticks([])

# Compute and plot second dendrogram, which will be beside the graph. This
# dendogram is for the metabolites.
ax2 = fig.add_axes([0,0.1,0.05,0.7], frame_on=False)
Y = sch.linkage(seedMatrix, method='average', metric='euclidean')
Z2 = sch.dendrogram(Y, orientation='right')
ax2.set_xticks([])
ax2.set_yticks([])


# Plot distance matrix.
axmatrix = fig.add_axes([0.05,0.1,0.8,0.8])
idx1 = Z1['leaves']
idx2 = Z2['leaves']
seedMatrix = seedMatrix[:,idx1]
seedMatrix = seedMatrix[idx2,:]
im = axmatrix.matshow(seedMatrix, aspect='auto', origin='lower')
axmatrix.set_xticks([])
axmatrix.set_yticks([])

axmatrix.set_xticks(range(len(seedMatrixT)))
axmatrix.set_xticklabels([dirList[i] for i in idx1], minor=False)
axmatrix.xaxis.set_label_position('bottom')
axmatrix.xaxis.tick_bottom()
pyplot.xticks(rotation=-90, fontsize=8)

axmatrix.set_yticks(range(len(seedMatrix)))
axmatrix.set_yticklabels([seedMatrixDF['Metabolite'][i] for i in idx2], minor=False)
axmatrix.yaxis.set_label_position('right')
axmatrix.yaxis.tick_right()

# Plot colorbar.
axcolor = fig.add_axes([0.97,0.1,0.03,0.8])
pyplot.colorbar(im, cax=axcolor)
fig.show()
fig.savefig('../'+summaryStatsDir+'/'+'seedSetDendrogram.png')
