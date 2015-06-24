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

# Preliminaries
# This cell imports the packags and modules used throughout this notebook and
# creates the list of genomes to process. The user is allowed to define their
# local folder structure for data input and processing.

# Import Python packages.
import cobra
import csv
import itertools
import os, glob
import networkx as nx
import pandas as pd

# Import scipy stack
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import scipy.cluster.hierarchy as sch
import graphFunctions as gf

# Define data directories
rawModelDir = 'ProcessedModelFiles'
processedDataDir = 'MergedData'
summaryStatsDir = 'DataSummaries/MergedData'
externalDataDir = 'ExternalData'

#%%

# The first step is to identify all samples belonging to the same tribe
# This function reads in a taxonomy and returns a dict containing the genomes
# associated with each tribe
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

#%%
    
# In the second step, all samples from a tribe are identified. Each sample is
    # converted to a graph object and merged with the previous graph. The final
    # graph is written to file.
    
def createTribalGraph(tribeSampleDict):

    print 'Merging genomes from individual tribes'
    
# Loop over the keys of the dictionary, one for each tribe
    for tribe in tribeSampleDict:
# Create an empty graph object
        tribalGraph = nx.DiGraph()
# Read in the graph of the tribe and merge with the graph from the previous
# iteration
        for sample in tribeSampleDict[tribe]:
# Read in adjacency list and convert to digraph object
            myDiGraph = nx.read_adjlist('../'+rawModelDir+'/'+sample+'/'+sample+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())
# Append to the previous graph
            tribalGraph = nx.compose(tribalGraph, myDiGraph)
# Fix edge weights to remove duplicate edges
# Find the largest SCC and reduce the graph to its SCC

# Check that the proper output directory exists. It not, create it.
        if not os.path.exists('../'+processedDataDir+'/'+tribe):
            os.makedirs('../'+processedDataDir+'/'+tribe)
    
        nx.write_adjlist(tribalGraph, '../'+processedDataDir+'/'+tribe+'/'+tribe+'AdjList.txt', delimiter='\t')

    return

#%%
def getDirList(inputDir):
    dirList =[]
    for item in os.listdir(inputDir):
        if not item.startswith('.'):
            dirList.append(item)
           
    return dirList
    
#%% Compute graph statistics. 
    
# This functions reads in the adjacency lists from the given directory and creates
# graph and directed graph (digraph) representations of each list. The objects
# are created using the networkX package. Summary statistics for the graph
# and directed graph are also reported and written to file.

def computeGraphStats(dirList):

# Create arrays to store summary statistics. Each array has four integer
# columns, for different properties of the graph: number of nodes (metabolites),
# edges, total components, and size of the largest component.
    
    graphStatArray = np.empty([numSubDir, 4], dtype = int)
    diGraphStatArray = np.empty([numSubDir, 4], dtype = int)

# Create files to record the summary statistics.
    graphFile = open('../'+summaryStatsDir+'/'+'GraphStatistics.txt', 'w')
    graphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

    diGraphFile = open('../'+summaryStatsDir+'/'+'DiGraphStatistics.txt', 'w')
    diGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

# Create lists to store the graph and digraph objects
    graphList = []
    diGraphList = []

# Iterate over the list of genome directories. For each genome, read in the
# adjacency list and convert it to both a graph and a digraph. Retrieve 
# properties of the (di)graph and record in the appropriate array. Write these
# statistics to file.
    count = 0
    print 'Computing Graph Statistics'

    for curDir in dirList:
# Read in adjacency list and convert to graph object
        myGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                              delimiter='\t', create_using=nx.Graph())

# Read in adjacency list and convert to digraph object
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())                            

# Append to the appropriate list
        graphList.append(myGraph)
        diGraphList.append(myDiGraph)
        
# Read model statistics by invoking graphFunctions.getGraphStats
        graphStatArray[count:] = gf.getGraphStats(myGraph)
        graphFile.write('%s,%i,%i,%i,%i\n' % (curDir, graphStatArray[count,0], 
                                       graphStatArray[count,1], 
                                       graphStatArray[count, 2],
                                       graphStatArray[count, 3] ) )

# Read model statistics by invoking graphFunctions.getDiGraphStats
        diGraphStatArray[count:] = gf.getDiGraphStats(myDiGraph)
        diGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, diGraphStatArray[count,0], 
                                       diGraphStatArray[count,1], 
                                       diGraphStatArray[count, 2],
                                       diGraphStatArray[count, 3] ) )

        count = count + 1
    
# Close files containing summary data
    graphFile.close()
    diGraphFile.close()
    
    return    
    
#%% 
# Reduction to largest component.

def reduceToLargeComponent(dirList):

# Create arrays to store summary statistics. Each array has four integer
# columns, for different properties of the graph: number of nodes (metabolites),
# edges, total components, and size of the largest component.
    reducedGraphStatArray = np.empty([numSubDir, 4], dtype = int)
    reducedDiGraphStatArray = np.empty([numSubDir, 4], dtype = int)

# Create files to record the summary statistics.
    reducedGraphFile = open('../'+summaryStatsDir+'/'+'ReducedGraphStatistics.txt', 'w')
    reducedGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

    reducedDiGraphFile = open('../'+summaryStatsDir+'/'+'ReducedDiGraphStatistics.txt', 'w')
    reducedDiGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

# Create lists to store the graph and digraph objects
    reducedGraphList = []
    reducedDiGraphList = []

# Iterate over the list of genome directories. For each graph, identify its
# largest component and identify the nodes belonging to each component.
# Aggregate the remaining nodes and remove them from both the graph and the
# digraph. Write the graph and diGraph to file as an adjacency list.
    count = 0
    print 'Reducing to Largest Component'

    for curDir in dirList:
    
# Read in adjacency list and convert to graph object
        myGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                              delimiter='\t', create_using=nx.Graph())

# Read in adjacency list and convert to digraph object
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())                            

# Identify the connected components of the graph representation and sort from
# largest to smallest (subGraphs). Aggregate the nodes in all but the largest 
# component (removeNodes).
        subGraphs = sorted(nx.connected_components(myGraph), key = len, reverse=True)
        removeNodes = list(itertools.chain(*subGraphs[1:len(subGraphs)]))

# Remove these nodes from both the graph and digraph    
        myGraph.remove_nodes_from(removeNodes)
        myDiGraph.remove_nodes_from(removeNodes)

# Append to the appropriate list
        reducedGraphList.append(myGraph)
        reducedDiGraphList.append(myDiGraph)

# Read model statistics by invoking graphFunctions.getGraphStats                
        reducedGraphStatArray[count:] = gf.getGraphStats(myGraph)
        reducedGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, reducedGraphStatArray[count,0], 
                                       reducedGraphStatArray[count,1], 
                                       reducedGraphStatArray[count, 2],
                                       reducedGraphStatArray[count, 3] ) )
# Read model statistics by invoking graphFunctions.getDiGraphStats
        reducedDiGraphStatArray[count:] = gf.getDiGraphStats(myDiGraph)
        reducedDiGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, reducedDiGraphStatArray[count,0], 
                                       reducedDiGraphStatArray[count,1], 
                                       reducedDiGraphStatArray[count, 2],
                                       reducedDiGraphStatArray[count, 3] ) )
# Create adjacency list for the reduced digraph and write to file
        nx.write_adjlist(myDiGraph, '../'+processedDataDir+'/'+curDir+'/'+curDir+'RedAdjList.txt')
                                       
        count = count + 1

# Close files containing summary data
    reducedGraphFile.close()
    reducedDiGraphFile.close()
    
    return
    
#%% 
    
def computeSeedSets(dirList):
        
# Computation of seed sets.

# This code cell does a number of things. First, it computes the strongly 
# connected components (SCCs) of the reduced digraph. An SCC is a group of
#  nodes, such that from each node there exists a path to all other nodes in 
# the component. SCCs are candidates for seed sets.

# Second, SCCs are evaluated to see if they are seed sets: any SCC with no 
# outgoing edges is a seed set. This is done by converting the digraph to its 
# condensation (a new graph) in which each SCC is represented as a single
# node. 

# Third, seed sets are written to file and summary statistics are computed
# for each seed set. Additional statistics on the reduced graph and digraph
# are also computed.

# Create lists to store seed sets
# seedSetList is a list of lists. Each outer list contains all the seed sets
# for that graph.
    seedSetList = []

# Iterate over the list of genome directories. For each reduced digraph, 
# identify its condensation (SCCs). For each node of the SCC, check if it
# is a seed set by computing its in-degree. If yes, append the SCC (as a list
# of nodes) to the list of seed sets. Then compute some summary statistics.
    count = 0
    print 'Computing Seed Sets'

    for curDir in dirList:

# Read in adjacency list and convert to digraph object
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())                            

# Compute the list of SCCs for the digraph as well as its condensation
        mySCCList = list(nx.strongly_connected_components_recursive(myDiGraph))
        myCondensation = nx.condensation(myDiGraph)

# "List of lists" of seed metabolites. Each element is a list of nodes belonging
# to an SCC which is also a seed set.
        mySeeds = []    

# For each node (SCC) of the condensation, examine each its in-degree. If the
# in-degree is zero (only outgoing edges), the SCC is a seed set. Append the
# SCC (as a list of nodes) to the list of seed sets.
        for node in myCondensation.nodes():
            inDeg = myCondensation.in_degree(node)
            if inDeg == 0:
                mySeeds.append(mySCCList[node])
            seedSetList.append(mySeeds)

# Record seed metabolites for each graph. Each row of the output file contains
# the metabolites belonging to a single seed set.
        seedSets = open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedSets.txt', 'w')
        writer = csv.writer(seedSets)
        writer.writerows(mySeeds)
        seedSets.close()
    
# Update the list of seed metabolites: replace the Model SEED metabolite 
# identifier with its common name. Note: The file metabMap.csv was created 
# manually from the seed database, and should be updated to reflect the 
# particulars of your data set. As above, record the seed metabolite for each
# graph.

# First read metabMap.csv in as a dictionary
        with open('../'+externalDataDir+'/'+'metabMap.csv', mode='rU') as inFile:
            reader = csv.reader(inFile)
            namesDict = dict((rows[0],rows[1]) for rows in reader)
        
# For each compound in the set of seeds, use the dictionary to replace it with
# its common name. Then write to file.
        mySeedsNames = [[namesDict[metab] for metab in seed] for seed in mySeeds]    

        seedSets = open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedSetsWNames.txt', 'w')
        writer = csv.writer(seedSets)
        writer.writerows(mySeedsNames)
        seedSets.close()
    
# Record weights for each seed metabolite. Each row of the output file contains
# a metabolite and its weight (1 / size of the seed set). Construct for seeds
# using both IDs and names.
        seedWeights = open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedWeights.txt', 'w')
        for seed in mySeeds:
            myWeight = 1 / float(len(seed))
            for metab in seed:
                seedWeights.write('%s,%f\n' % (metab, myWeight) )
        seedWeights.close()
    
        seedWeights = open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedWeightsWNames.txt', 'w')
        for seed in mySeedsNames:
            myWeight = 1 / float(len(seed))
            for metab in seed:
                seedWeights.write('%s,%f\n' % (metab, myWeight) )
        seedWeights.close()
    
        count = count + 1

    return


#%% Consolidation seed weights into a single data frame.
# This code snippet reads in the seed metabolites and their weights for each 
# genome. The data is read into a dataframe and then written to file. Structure:
# Rows: metabolites
# Columns: graphs
# Entries: unweighted seed set values

def consolidateSeeds(dirList):
    
    print 'Consolidate seed sets'
        
# Create a data frame to store the data
    seedMatrixDF = pd.DataFrame(columns=["Metabolite"])

# Read in list of seed weights into a temporary data frame. Perform an outer
# join with the existing data frame to incorporate the new list of weights.
    for curDir in dirList:
        tempDF = pd.read_csv('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedWeights.txt', names=['Metabolite',curDir])
        seedMatrixDF = pd.merge(seedMatrixDF, tempDF, how='outer', on="Metabolite")

# Replace all the NaN values with zeros
    seedMatrixDF.fillna(0, inplace=True)

# Append a new column containing common names associated with metabolite IDs.
# The file metabMap.csv was created manually from the seed database, and should
# be updated to reflect the particulars of your data set.
    namesDF = pd.read_csv('../'+externalDataDir+'/'+'metabMap.csv', names=['Metabolite','CommonName'])
    seedMatrixDF = pd.merge(seedMatrixDF, namesDF, how='inner', on="Metabolite")

# Rearrange the order of the columns so that the common name is in front
    newOrder = seedMatrixDF.columns.tolist()
    newOrder = newOrder[-1:] + newOrder[:-1]
    seedMatrixDF = seedMatrixDF[newOrder]
    
# Export the matrix of seed weights
    pd.DataFrame.to_csv(seedMatrixDF, '../'+summaryStatsDir+'/'+'seedMatrixWeighted.csv')
    
    return seedMatrixDF


#%% Clustering and Visualization of Seed Sets

def clusterSeedSets(seedMatrixDF):
# This code snippet creates a dendrogram of the seed sets. The genomes are
# clustered based on weighted vectors of seed sets, using the euclidean
# distance and UPGMA (average linkage) clustering. A second dendrogram is
# constructed based on metabolite weights across genomes. Then, the matrix
# of seed weights is reordered to reflect the order of the dendrograms and is
# visualized. Genome names are colored according to their lineage, with acI
# sub-divided into acI-A and acI-B.

    print 'Computing dendrogram'
    
# Python clustering algorithms require the data to be an ndarray, with each
# row corresponding to a set of observations. Thus, calculation of the linkage
# for the genomes must be performed on the transpose.
    seedMatrix=pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1))
    seedMatrixT=np.transpose(pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1)))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7.5, 61))

# Compute and plot dendrogram for genomes, which will be above the graph.
# Define the size of the dendrogram
    ax1 = fig.add_axes([0.05,0.9,0.8,0.05], frame_on=False)
# Compute the linkage matrix
    genomeLinkage = sch.linkage(seedMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    genomeClust = sch.dendrogram(genomeLinkage, distance_sort='True', color_threshold=0)
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Compute and plot dendrogram for metabolites which will be beside the graph. 
# Define the size of the dendrogram
    ax2 = fig.add_axes([0,0.1,0.05,0.8], frame_on=False)
# Compute the linkage matrix
    metabLinkage = sch.linkage(seedMatrix, method='average', metric='euclidean')
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
    seedMatrix = seedMatrix[:,idx1]
    seedMatrix = seedMatrix[idx2,:]
# Plot the weight matrix
    im = axmatrix.matshow(seedMatrix, aspect='auto', origin='lower')
# No tick marks along axes
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+'tribalColors.csv')
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]

# Add genome names to the bottom axis
    axmatrix.set_xticks(range(len(seedMatrixT)))
    axmatrix.set_xticklabels([dirList[i] for i in idx1], minor=False)
    axmatrix.xaxis.set_label_position('bottom')
    axmatrix.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(axmatrix.get_xticklabels(), genomeColors):
        xtick.set_color(color)

# Add metabolite names to the right axis
    axmatrix.set_yticks(range(len(seedMatrix)))
    axmatrix.set_yticklabels([seedMatrixDF['CommonName'][i] for i in idx2], minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    plt.yticks(fontsize=8)

# Plot colorbar.
    axcolor = fig.add_axes([0.97,0.9,0.03,0.05])
    plt.colorbar(im, cax=axcolor)
    fig.show()
    fig.savefig('../'+summaryStatsDir+'/'+'seedSetDendrogram.png')

#%% Clustering and Visualization of Seed Sets

def clusterOnly(seedMatrixDF):

    seedMatrix=pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1))
    seedMatrixT=np.transpose(pd.DataFrame.as_matrix(seedMatrixDF.drop(['CommonName','Metabolite'], 1)))

# Create a figure to display the seed weights and dendrograms.
    fig = plt.figure(figsize=(7,2))
    
# Define the size of the dendrogram
    ax1 = fig.add_axes([0, 0, 1, 1], frame_on=False)
# Compute the linkage matrix
    genomeLinkage = sch.linkage(seedMatrixT, method='average', metric='euclidean')
# Compute the dendrogram
    genomeClust = sch.dendrogram(genomeLinkage, distance_sort='True', color_threshold=0)
    idx1 = genomeClust['leaves']
# No tick marks along axes
    ax1.set_xticks([])
    ax1.set_yticks([])

# Import coloration info to map to genome names. Rearrange to same order as
# leaves of the dendrogram and extract the 'Color' column as a list.
# The file 'actinoColors.csv' will need to be updated for the specific samples.
    genomeColors = pd.read_csv('../'+externalDataDir+'/'+'tribalColors.csv')
    genomeColors = genomeColors['Color'].tolist()
    genomeColors = [ genomeColors[i] for i in idx1]
    
    axisLength = 1- (1/float(len(dirList)))
    axisStart = (1 - axisLength) / 2
# Add genome names to the bottom axis
    ax2 = fig.add_axes([axisStart, 0, axisLength, 0], frame_on=False)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xticks(range(len(seedMatrixT)))
    ax2.set_xticklabels([dirList[i] for i in idx1], minor=False)
    ax2.xaxis.set_label_position('bottom')
    ax2.xaxis.tick_bottom()
    plt.xticks(rotation=-90, fontsize=8)
    for xtick, color in zip(ax2.get_xticklabels(), genomeColors):
        xtick.set_color(color)


#%% Actual program operation
# Will go at end of file. Temporarily moved above working code to ensure all
# prereq functions have been called.    
tribeSampleDict =  importTaxonomy('../ExternalData/taxonomySAGs.csv')
createTribalGraph(tribeSampleDict)
dirList = getDirList('../'+processedDataDir)
numSubDir = len(dirList)
computeGraphStats(dirList)
reduceToLargeComponent(dirList)
computeSeedSets(dirList)
seedMatrixDF = consolidateSeeds(dirList)
clusterSeedSets(seedMatrixDF)
clusterOnly(seedMatrixDF)