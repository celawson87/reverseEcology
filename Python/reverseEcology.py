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

# Compute statistics on the size of metabolic network graphs
gf.computeGraphStats(dirList, processedDataDir, summaryStatsDir)

# Reduce network graphs to their largest component
gf.reduceToLargeComponent(dirList, processedDataDir, summaryStatsDir)

# Compute seed sets
gf.computeSeedSets(dirList, externalDataDir, processedDataDir)

# Aggregate seed sets for all models and write to file
seedMatrixDF = ef.consolidateSeeds(dirList, externalDataDir, processedDataDir, summaryStatsDir)

# Calculate normalized counts of seed metabolites
ef.normalizedSeedCounts(dirList, processedDataDir, summaryStatsDir)

# Construct heatmap of seed metabs. Cluster seed profiles and create dendrogram
ef.clusterSeedSets(seedMatrixDF, dirList, externalDataDir, summaryStatsDir, 'actinoColors.csv')
ef.clusterOnly(seedMatrixDF, dirList, externalDataDir, 'actinoColors.csv')

# Compute metabolic competition and cooperation scores. Create heatmap. Cluster 
# each metric and create a dendrogram
metabCompeteDF = ef.computeMetabCompete(dirList, processedDataDir, summaryStatsDir)
ef.clusterPairwise(metabCompeteDF, dirList, externalDataDir, summaryStatsDir, 'actinoColors.csv', 'metabolicCompetition.png')

metabComplementDF = ef.computeMetabComplement(dirList, processedDataDir, summaryStatsDir)
ef.clusterPairwise(metabComplementDF, dirList, externalDataDir, summaryStatsDir, 'actinoColors.csv', 'metabolicComplementarity.png')
