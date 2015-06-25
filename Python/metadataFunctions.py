###############################################################################
# metadataFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with metadata and other external datasets.
################################################################################

# Import Python packages.
import os, glob
import pandas as pd

################################################################################

# getDirList
# Retrieve list of genomes to process by examing the contents of 'inputDir,'
# ignoring hidden folders. Each folder contains data files for that genome 
# (e.g., metabolic models, network graphs, etc), and calculations are performed
# performed by iterating over this list.

def getDirList(inputDir):
    dirList =[]
    for item in os.listdir(inputDir):
        if not item.startswith('.'):
            dirList.append(item)
            
    return dirList
    
################################################################################

# importTaxonomy
# This function reads in a taxonomy file and returns a dicttionary containing 
# the genomes associated with each tribe. The file should a csv file organized
# as follows:
# Sample	Lineage	Clade	Tribe
# AAA023D18	acI	acI-B	acI-B1
# The function could be updated to take (lineage, clade, tribe) as input

def importTaxonomy(taxonFile):
    
    print 'Importing taxonomy'

# Read in the taxonomic classification
    taxonClass = pd.DataFrame.from_csv(taxonFile, sep=',')
    taxonClass = taxonClass.dropna()
    
# Extract the unique tribes found in the dataset
    tribes = pd.unique(taxonClass.Tribe.values)
    tribes.sort(axis=0)
    
# For each tribe, return the list of samples. Creates a dict and adds an entry
    # for each tribe.
    tribeSampleDict = {}

    for tribe in tribes:
# Identify the samples belonging to this tribe
        samples = taxonClass.loc[taxonClass['Tribe'] == tribe]
        samples = [sample for sample in samples.index]
        tribeSampleDict[tribe] = samples
        
    return tribeSampleDict