###############################################################################
# markerGeneFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for extracting marker gene presence/absence from Phylosift
# and Phylophlan.
################################################################################ 

# Import Python packages.
import os, glob
import pandas as pd
import re

################################################################################

# Extract marker gene presence/absence from phylosift data. Use specifies the
# path to phylosift and an output file name. Matrix of presence/absence of
# marker genes is written as a CSV file.
def extractPhylosiftMarkerGenes(inputDir, fileName):
# List all the folders in the directory
    dirList =[]
    for item in os.listdir(inputDir+'/PS_temp'):
        if not item.startswith('.'):
            dirList.append(item)
        
# Initialize an empty dataframe
        markerGeneDF = pd.DataFrame(columns=["MarkerGene"])

# Read in the markerSummary file as a dataframe and append to the previous        
    for curDir in dirList:
        tempDF = pd.read_csv(inputDir+'/PS_temp/'+curDir+'/alignDir/marker_summary.txt', names=['MarkerGene',curDir], delimiter='\t')
        markerGeneDF = pd.merge(markerGeneDF, tempDF, how='outer', on="MarkerGene")

# Clean up the final data frame.
# Drop unwanted marker genes
    markerGeneDF = markerGeneDF[markerGeneDF.MarkerGene.str.contains('DNGNGWU')]
# Sort alphabetically
    markerGeneDF = markerGeneDF.sort(['MarkerGene'])
# Replace all the NaN values with zeros
    markerGeneDF.fillna(0, inplace=True)
# Fix column names, removing '.fna', 'metabat', and 'v2'
    markerGeneDF = markerGeneDF.rename(columns=lambda x: re.sub('.fna','',x))
    markerGeneDF = markerGeneDF.rename(columns=lambda x: re.sub('metabat','',x))
    markerGeneDF = markerGeneDF.rename(columns=lambda x: re.sub('v2','',x))

# Write to file, without the index column
    pd.DataFrame.to_csv(markerGeneDF, fileName, index=False)

    return

################################################################################
# Extract marker gene presence/absence from phylosift data. Use specifies the
# path to phylosift and an output file name. Matrix of presence/absence of
# marker genes is written as a CSV file.

def extractPhylophlanMarkerGenes(inputDir, fileName):
# Grab a list of *.b6o files which contain marker gene info
    fileList = [os.path.basename(x) for x in glob.glob(inputDir+'/*.b6o')]

# Initialize an empty dataframe
    markerGeneDF = pd.DataFrame(columns=["MarkerGene"])

# Read in the markerSummary file as a dataframe and append to the previous        
    for curFile in fileList:
        tempDF = pd.read_csv(inputDir+'/'+curFile, names=['Locus','MarkerGene'], delimiter='\t')
# Drop the first column and append a second column with all '1s' and genome name as the header
        tempDF = tempDF.drop('Locus', axis=1)
        tempDF[curFile] = 1

        markerGeneDF = pd.merge(markerGeneDF, tempDF, how='outer', on="MarkerGene")
    
# Clean up the final data frame.
# Sort alphabetically
        markerGeneDF = markerGeneDF.sort(['MarkerGene'])
# Replace all the NaN values with zeros
        markerGeneDF.fillna(0, inplace=True)
# Fix column names, removing '.fna', 'metabat', and 'v2'
        markerGeneDF = markerGeneDF.rename(columns=lambda x: re.sub('.b6o','',x))
        markerGeneDF = markerGeneDF.rename(columns=lambda x: re.sub('metabat','',x))
        markerGeneDF = markerGeneDF.rename(columns=lambda x: re.sub('v2','',x))

# Write to file, without the index column
    pd.DataFrame.to_csv(markerGeneDF, fileName, index=False)

    return

################################################################################
# Function Calls
extractPhylosiftMarkerGenes('phylosift', 'phylosiftMarkerGenes.csv')
extractPhylophlanMarkerGenes('phylophlan/actino', 'phylophlanMarkerGenes.csv')
