###############################################################################
# pairwiseANI.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Input:
#   aniFile - table of pairwise ANI or coverage comparisons
#   taxon.csv - table of samples w/ taxonomic classification: lineage, clade, tribe
# Output: 
#   Dataframe giving the following information:
#   tribe   samples   min pairwiseANI/COV   max pairwiseANI/COV
################################################################################

# Import python packages
import itertools
    
# Import scipy stack
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy

#%% 

def importANIandTaxonomy(externalDataDir, aniFile, taxonFile):
# Read in the pairwise ANI calculations. Replace '-' with NaN so we can tell
# pandas to ignore it.
    pairwiseANI = pd.DataFrame.from_csv('../'+externalDataDir+'/'+aniFile, sep='\t')
    pairwiseANI = pairwiseANI.convert_objects(convert_numeric=True)

# Read in the taxonomic classification
    taxonClass = pd.DataFrame.from_csv('../'+externalDataDir+'/'+taxonFile, sep=',')
    taxonClass = taxonClass.dropna()

# Update the pairwiseANI dataframe to only include samples with a taxonomic
# classification to the tribe level
    samples = taxonClass.index
    redPairwiseANI = pairwiseANI.loc[samples,samples]
    redPairwiseANI.to_csv('../'+externalDataDir+'/'+'reduced'+aniFile)

# Extract the unique tribes found in the dataset
    tribes = pd.unique(taxonClass.Tribe.values)
    tribes.sort(axis=0)

    return pairwiseANI, taxonClass, tribes

#%% Calculate min and max pairwise ANI for members of the same tribe

def sameTribePairwiseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile):
# Create a dataframe to store output. Indices are tribes, columns are the
# the data assoc. w/ each tribe
    maxMinANI = pd.DataFrame(index = tribes, columns=['Samples', 'Num Samples', 'Max', 'Min'])

# Iterate over the set of unique tribes. Create a smaller dataFrame containing
# just the pairwiseANI scores from that tribe. Return the max and min and
# insert them into the dataframe. When complete, write the dataframe to file.

    for tribe in tribes:
# Identify the samples belonging to this tribe
        samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
        samples = [sample for sample in samples.index]

# Extract the reduced dataframe and write to file
        redPairwiseANI = pairwiseANI.loc[samples,samples]
#    redPairwiseANI.to_csv(tribe+'.csv')

# Compute the max and min ANI
        maxANI = redPairwiseANI.max().max()
        minANI = redPairwiseANI.min(skipna=True).min()
    
# Add this information to the dataframe
        maxMinANI['Samples'][tribe] = samples
        maxMinANI['Num Samples'][tribe] = len(samples)
        maxMinANI['Max'][tribe] = maxANI
        maxMinANI['Min'][tribe] = minANI
        maxMinANI.to_csv('../'+externalDataDir+'/'+'withinTribe'+aniFile+'.csv')
        
    return maxMinANI
    
#%% Calculate min and max pairwise ANI for members of differing tribes

def diffTribePairwiseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile):
# Create a dataframe to store output. Indices are tribes, columns are the
# the data assoc. w/ each tribe
    maxMinANI = pd.DataFrame(index = tribes, columns=['Max', 'Min'])

# Iterate over the set of unique tribes. Identify samples belonging to that
# tribe and remainint tribes. Create a smaller dataFrame containing just the 
# pairwiseANI scores from those samples. Return the max and min and
# insert them into the dataframe. When complete, write the dataframe to file.
    for innerTribe in tribes:
        innerSamples = taxonClass.loc[taxonClass['Tribe'] == innerTribe]
        innerSamples = [innerSample for innerSample in innerSamples.index]
        outerSamples = taxonClass.loc[taxonClass['Tribe'] != innerTribe]
        outerSamples = [outerSample for outerSample in outerSamples.index]
    
# Extract the reduced dataframe and write to file
        redPairwiseANI = pairwiseANI.loc[innerSamples,outerSamples]
#        redPairwiseANI.to_csv(innerTribe+'-'+outerTribe+'.csv')

# Compute the max and min ANI
        maxANI = redPairwiseANI.max().max()
        minANI = redPairwiseANI.min(skipna=True).min()

# Add this information to the dataframe
        maxMinANI['Max'][innerTribe] = maxANI
        maxMinANI['Min'][innerTribe] = minANI
        maxMinANI.to_csv('../'+externalDataDir+'/'+'betweenTribe'+aniFile+'.csv')
        
    return maxMinANI
    
#%% Determine if a new genome can be added to a group of samples belonging
# to a single tribe
    
def addGenomeToTribe(pairwiseANI, taxonClass, tribes, tribe, newGenomes):
# Compute max and min pairwise ANI for all samples in the same tribe
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]
    redPairwiseANI = pairwiseANI.loc[samples,samples]
    minANItribe = redPairwiseANI.min(skipna=True).min()

# Now determine the pairwise ANI between the new sample and the tribe of
# interest. Append the new sample to the sample list and update the set of
# pairwise ANI.
    samples.append(newGenomes)
    redPairwiseANI = pairwiseANI.loc[samples,samples]
    minANI = redPairwiseANI.min(skipna=True).min()

    print redPairwiseANI
    print str(minANItribe)

    print('\n'+'The minimum within the tribe '+tribe+' is: '+str(minANItribe))
    print('\n'+'When genome '+newGenomes+' is added, the new minimum becomes: '+str(minANI))

    return

#%% Determine if a new genome can be added to a single sample based on the
# smallest pairwise ANI from any two samples belonging to the same trobe

def compareSamples(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile, existingGenomes, newGenome):
# Compute max and min pairwise ANI for all samples in the same tribe
    maxMinANI = sameTribePairwiseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile)

# Compute the overall minANI
    universalMin = maxMinANI['Min'].min()
        
# Compute maximal pairwise ANI of sample being considered for merger
    samples = existingGenomes + newGenome
    redPairwiseANI = pairwiseANI.loc[samples,samples]
    minANI = redPairwiseANI.min(skipna=True).min()

#    print('\n'+'The smallest ANI cutoff for any tribe is: '+str(universalMin))
    print('\n'+'When genome '+str(newGenome)+' is added, the min among all samples: '+str(minANI))
    
    return
    
#%% For a given tribe and sample size, determine all possible combinations
# of SAGs of that sample size. For each, compute their maximum and minimum 
# pairwise ANI. Report the worst-case scenario, i.e., the highest minimum.
    
def worstCaseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile, tribe, sampleSize):
    
    sameTribePairwiseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile)

# Create the list of samples
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]

# Create a list of all combinations of specified size
    randomSamples = list(itertools.combinations(samples, sampleSize))

# Create a dataframe to store the results
    maxMinANI = pd.DataFrame(index = randomSamples, columns=['Max', 'Min'])

# Loop over all samples and calculate min and max pairwise ANI
    for sample in randomSamples:            
        redPairwiseANI = pairwiseANI.loc[sample,sample]

# Compute the max and min ANI
        maxANI = redPairwiseANI.max().max()
        minANI = redPairwiseANI.min(skipna=True).min()

# Add this information to the dataframe
        maxMinANI['Max ANI'][sample] = maxANI
        maxMinANI['Min ANI'][sample] = minANI

    print maxMinANI
    print('Sampling tribe '+tribe+ ' with sample size ' +str(sampleSize)+ '.')
    print('The worst-case minimum ANI is: ' +str(maxMinANI['Min'].max()))

    return
    
#%% For a given tribe, compute the worst case scenario ANI for each
# possible sample size.

def allWorstCaseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile, tribe):
    
    sameTribePairwiseANI(externalDataDir, pairwiseANI, taxonClass, tribes, aniFile)

# Create the list of samples
    samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
    samples = [sample for sample in samples.index]

# Create a numpy array to store the worst-case value for each pair
    sampleMin = np.empty([2, len(samples)-1,])

    for count in range (2, len(samples)+1):
# Create a list of all combinations of specified size
        randomSamples = list(itertools.combinations(samples, count))

# Create a dataframe to store the results
        maxMinANI = pd.DataFrame(index = randomSamples, columns=['Max', 'Min'])

# Loop over all samples and calculate min and max pairwise ANI
        for sample in randomSamples:            
            redPairwiseANI = pairwiseANI.loc[sample,sample]

# Compute the max and min ANI
            maxANI = redPairwiseANI.max().max()
            minANI = redPairwiseANI.min(skipna=True).min()

# Add this information to the dataframe
            maxMinANI['Max'][sample] = maxANI
            maxMinANI['Min'][sample] = minANI

        sampleMin[0][count-2] = count
        sampleMin[1][count-2] = maxMinANI['Min'].max()
        print('Sampling tribe '+tribe+ ' with sample size ' +str(count)+ '.')
        print('The worst-case minimum is: ' +str(maxMinANI['Min'].max()))
    
    plt.scatter(sampleMin[0], sampleMin[1])
    plt.xlim(1, len(samples)+1)
    plt.ylim(0.9*sampleMin[1].min(), 100)
    plt.xlabel('Number of Samples')
    plt.ylabel('Cutoff for Tribe')
    
    return
