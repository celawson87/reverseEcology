###############################################################################
# meringGraphs.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# This notebook creates "merged" metabolic reconstructions by aggregating
# graphs for all SAGs from the same tribe,
################################################################################

# Import Python modules 
# These custom-written modules should have been included with the package
# distribution. 
import graphFunctions as gf
import metadataFunctions as mf
import sbmlFunctions as sf
import seedFunctions as ef

# Define local folder structure for data input and processing.
rawModelDir = 'ProcessedModelFiles'
processedDataDir = 'MergedData'
summaryStatsDir = 'DataSummaries/MergedData'
externalDataDir = 'ExternalData'
        
#%% Program execution

# Identify all samples belonging to the same tribe
tribeSampleDict =  mf.importTaxonomy('../ExternalData/taxonomySAGs.csv')

# Create a graph for each tribe by merging graphs from individual members
# of that tribe
gf.createTribalGraph(tribeSampleDict, processedDataDir, rawModelDir)

# Import the list of models
dirList = mf.getDirList('../'+processedDataDir)
numSubDir = len(dirList)

# Compute graph statistics
gf.computeGraphStats(dirList, processedDataDir, summaryStatsDir)

# Reduction to largest component
gf.reduceToLargeComponent(dirList, processedDataDir, summaryStatsDir)

# Compute seed sets
gf.computeSeedSets(dirList, externalDataDir, processedDataDir)

# Aggregate seed sets for all models and write to file
seedMatrixDF = ef.consolidateSeeds(dirList, externalDataDir, processedDataDir, summaryStatsDir)

# Construct heatmap of seed metabs. Cluster seed profiles and create dendrogram
ef.clusterSeedSets(seedMatrixDF, dirList, externalDataDir, summaryStatsDir)
ef.clusterOnly(seedMatrixDF, dirList, externalDataDir)
