###############################################################################
# reverseEcology.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Master function for performing reverse ecology simulations.
################################################################################

# Import Python modules 
# These custom-written modules should have been included with the package
# distribution. 
import graphFunctions as gf
import metadataFunctions as mf
import sbmlFunctions as sf
import seedFunctions as ef

# Define local folder structure for data input and processing.
processedDataDir = 'ProcessedModelFiles'
summaryStatsDir = 'DataSummaries'
externalDataDir = 'ExternalData'

#%% Program execution

# Import the list of models
dirList = mf.getDirList('../'+processedDataDir)
numSubDir = len(dirList)

# Convert SBML model files to adjacency lists
sf.dirListToAdjacencyList(dirList, externalDataDir, processedDataDir, summaryStatsDir)

# Compute graph statistics
gf.computeGraphStats(dirList, processedDataDir, summaryStatsDir)

# Reduction to largest component
gf.reduceToLargeComponent(dirList, processedDataDir, summaryStatsDir)

# Compute seed sets
gf.computeSeedSets(dirList, externalDataDir, processedDataDir)