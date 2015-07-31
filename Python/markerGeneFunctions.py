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

################################################################################
    
def importCogAndTaxonomy(externalDataDir, cogFile, taxonFile):
# Read in the COG abundances.
    cogDF = pd.DataFrame.from_csv('../'+externalDataDir+'/'+cogFile, sep=',')

# Read in the taxonomic classification
    taxonClass = pd.DataFrame.from_csv('../'+externalDataDir+'/'+taxonFile, sep=',')
    taxonClass = taxonClass.dropna()

# Extract the unique tribes found in the dataset
    tribes = pd.unique(taxonClass.Tribe.values)
    tribes.sort(axis=0)

    return cogDF, taxonClass, tribes

################################################################################    

def worstCasePanGenome(externalDataDir, cogAbundance, taxonClass, tribes, tribe, sampleSize):

# Create the list of samples
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]

# Create a list of all combinations of specified size
    randomSamples = list(itertools.combinations(samples, sampleSize))
    
# Create a dataframe to store the results
    corePanSize = pd.DataFrame(index = randomSamples, columns=['Core Genome', 'Pan Genome'])

# Loop over all samples and calculate core- and pan-genome
# Core-genome will be computed using an inner join
# Pan-genome will be computed using an outer join
# Because samples can be of arbitrary size, need to use a loop

    for sample in randomSamples:            
# Create the reduced dataframe
        redCogAbundance = cogAbundance[list(sample)]
    
# Use the value_counts function to count the number of instances where a COG
# is present in ALL genomes (core genome) or ANY genome (pangenome)
        coreGenome = (redCogAbundance != 0).all(axis=1).value_counts()
        panGenome = (redCogAbundance != 0).any(axis=1).value_counts()

# Add this information to the dataframe
        corePanSize['Core Genome'][sample] = coreGenome.loc[True]
        corePanSize['Pan Genome'][sample] = panGenome.loc[True]

    print corePanSize
    print('Sampling tribe '+tribe+ ' with sample size ' +str(sampleSize)+ '.')
    print('The smallest pan genome is: ' +str(corePanSize['Pan Genome'].min()))
    print('The largest pan genome is: ' +str(corePanSize['Pan Genome'].max()))

    return
    
################################################################################
    
def allWorstCasePanGenome(externalDataDir, cogAbundance, taxonClass, tribes, tribe):

# Create the list of samples
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]

# Create a numpy array to store the worst-case value for each sample size
    sampleMin = np.empty([3, len(samples)-1,])

# Loop over the range of possible sample sizes
    for sampleSize in range (2, len(samples)+1):
        
# Create a list of all combinations of specified size
        randomSamples = list(itertools.combinations(samples, sampleSize))
    
# Create a dataframe to store the results
        corePanSize = pd.DataFrame(index = randomSamples, columns=['Count', 'Core Genome', 'Pan Genome'])

# Loop over all samples and calculate core- and pan-genome
# Core-genome will be computed using an inner join
# Pan-genome will be computed using an outer join
# Because samples can be of arbitrary size, need to use a loop

        for sample in randomSamples:            
# Create the reduced dataframe
            redCogAbundance = cogAbundance[list(sample)]

# Use the value_counts function to count the number of instances where a COG
# is present in ALL genomes (core genome) or ANY genome (pangenome)
            coreGenome = (redCogAbundance != 0).all(axis=1).value_counts()
            panGenome = (redCogAbundance != 0).any(axis=1).value_counts()

# Add this information to the dataframe
            corePanSize['Count'][sample] = sampleSize
            corePanSize['Core Genome'][sample] = coreGenome.loc[True]
            corePanSize['Pan Genome'][sample] = panGenome.loc[True]
    
# Add this info to the master dataframe
        sampleMin[0][sampleSize-2] = sampleSize
# Find the min pangenome size
        sampleMin[2][sampleSize-2] = corePanSize['Pan Genome'].min()
# Find the index associated with this value and store the corresponding core
# genome size
        index = corePanSize['Pan Genome'].idxmin()
        sampleMin[1][sampleSize-2] = corePanSize['Core Genome'][index]
    
# Plot the results    
        cg = plt.scatter(sampleMin[0], sampleMin[1], color = 'b')
        pg = plt.scatter(sampleMin[0], sampleMin[2], color = 'g')
        plt.legend((cg, pg), ('Core Genome', 'Pan Genome'), loc='upper left')
        plt.xlim(1, len(samples)+1)
        plt.ylim(0, 1.1*sampleMin[2].max())
        plt.xlabel('Number of Samples')
        plt.ylabel('Core- or Pangenome Size')
    return