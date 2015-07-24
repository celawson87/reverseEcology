###############################################################################
# markerGeneFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with marker genes.
################################################################################ 

# Import python packages
import itertools
#from collections import Counter
#import csv
#import networkx as nx

# Import scipy stack
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#import scipy.cluster.hierarchy as sch

################################################################################

def importMarkerGenesAndTaxonomy(externalDataDir, markerGeneFile, taxonFile):
# Read in the pairwise ANI calculations. Replace '-' with NaN so we can tell
# pandas to ignore it.
    markerGeneDF = pd.DataFrame.from_csv('../'+externalDataDir+'/'+markerGeneFile, sep=',')
#    markerGeneDF = pairwiseANI.convert_objects(convert_numeric=True)

# Read in the taxonomic classification
    taxonClass = pd.DataFrame.from_csv('../'+externalDataDir+'/'+taxonFile, sep=',')
    taxonClass = taxonClass.dropna()

# Extract the unique tribes found in the dataset
    tribes = pd.unique(taxonClass.Tribe.values)
    tribes.sort(axis=0)

    return markerGeneDF, taxonClass, tribes
    
################################################################################

def indivCompleteness(externalDataDir, markerGeneDF, totalGenes, taxonClass, tribes):
# Create a dataframe to store output. Indices are tribes, columns are the
# the data assoc. w/ each tribe
    maxMinCompleteness = pd.DataFrame(index = tribes, columns=['Samples', 'Num Samples', 'Max', 'Min'])

# Iterate over the set of unique tribes. Create a smaller dataFrame containing
# just the completeness scores from that tribe. Return the max and min and
# insert them into the dataframe. When complete, write the dataframe to file.

    for tribe in tribes:
# Identify the samples belonging to this tribe
        samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
        samples = [sample for sample in samples.index]

# Extract the reduced dataframe and write to file
        redMarkerGeneDF = markerGeneDF[samples]
        indivCompleteness = redMarkerGeneDF.sum(axis = 0) / totalGenes
# Compute the completeness of each tribe

# Compute the max and min ANI
        maxCompleteness = indivCompleteness.max().max()
        minCompleteness = indivCompleteness.min(skipna=True).min()
    
# Add this information to the dataframe
        maxMinCompleteness['Samples'][tribe] = samples
        maxMinCompleteness['Num Samples'][tribe] = len(samples)
        maxMinCompleteness['Max'][tribe] = maxCompleteness
        maxMinCompleteness['Min'][tribe] = minCompleteness
        
    return maxMinCompleteness
        
################################################################################

def worstCaseCompleteness(externalDataDir, markerGeneDF, totalGenes, taxonClass, tribes, tribe, sampleSize):
    
    indivCompleteness(externalDataDir, markerGeneDF, totalGenes, taxonClass, tribes)

# Create the list of samples
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]

# Create a list of all combinations of specified size
    randomSamples = list(itertools.combinations(samples, sampleSize))
    
# Create a dataframe to store the results
    completeness = pd.DataFrame(index = randomSamples, columns=['Count', 'Percent'])

# Loop over all samples and calculate min and max pairwise ANI
    for sample in randomSamples:            
        redMarkerGeneDF = markerGeneDF[list(sample)]
        mergedDF = redMarkerGeneDF.sum(axis = 1)

# Compute the max and min ANI
        count = sum(mergedDF != 0)  
        percent = count / float(totalGenes)

# Add this information to the dataframe
        completeness['Count'][sample] = count
        completeness['Percent'][sample] = percent

    print completeness
    print('Sampling tribe '+tribe+ ' with sample size ' +str(sampleSize)+ '.')
    print('The worst-case estimated completeness is: ' +str(completeness['Percent'].min()))

    return
    

################################################################################

def allWorstCaseCompleteness(externalDataDir, markerGeneDF, totalGenes, taxonClass, tribes, tribe):
    
    indivCompleteness(externalDataDir, markerGeneDF, totalGenes, taxonClass, tribes)

# Create the list of samples
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]

# Create a numpy array to store the worst-case value for each pair
    sampleMin = np.empty([2, len(samples)-1,])

    for sampleSize in range (2, len(samples)+1):

# Create a list of all combinations of specified size
        randomSamples = list(itertools.combinations(samples, sampleSize))
    
# Create a dataframe to store the results
        completeness = pd.DataFrame(index = randomSamples, columns=['Count', 'Percent'])

# Loop over all samples and calculate min and max pairwise ANI
        for sample in randomSamples:            
            redMarkerGeneDF = markerGeneDF[list(sample)]
            mergedDF = redMarkerGeneDF.sum(axis = 1)

# Compute the max and min ANI
            count = sum(mergedDF != 0)  
            percent = count / float(totalGenes)

# Add this information to the dataframe
            completeness['Count'][sample] = count
            completeness['Percent'][sample] = percent

        sampleMin[0][sampleSize-2] = sampleSize
        sampleMin[1][sampleSize-2] = completeness['Percent'].min()
        print('Sampling tribe '+tribe+ ' with sample size ' +str(sampleSize)+ '.')
        print('The worst-case completeness is: ' +str(completeness['Percent'].min()))
    
    plt.scatter(sampleMin[0], sampleMin[1])
    plt.xlim(1, len(samples)+1)
    plt.ylim(0.9*sampleMin[1].min(), 1)
    plt.xlabel('Number of Samples')
    plt.ylabel('EstimatedCompleteness')

    return