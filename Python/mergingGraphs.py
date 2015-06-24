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

#%% Compute graph statistics. 
    
# This functions reads in the adjacency lists from the given directory and creates
# graph and directed graph (digraph) representations of each list. The objects
# are created using the networkX package. Summary statistics for the graph
# and directed graph are also reported and written to file.

def computeGraphStats(inputDir):

    dirList =[]
    for item in os.listdir(inputDir):
        if not item.startswith('.'):
            dirList.append(item)
            
    numSubDir = len(dirList)

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
        print 'Processing directory', count+1, 'of', numSubDir, ':', curDir

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
    
# Plot summary statistics. The function plotGraphStats plots histograms of:
#   graph size (number of nodes)
#   total number of components
#   size of largest compmonent, as fraction of total nodes
    gf.plotGraphStats(graphStatArray)
    
    return    
    
#%% 
# Reduction to largest component.

def reduceToLargeComponent(inputDir):
    
    dirList =[]
    for item in os.listdir(inputDir):
        if not item.startswith('.'):
            dirList.append(item)
            
    numSubDir = len(dirList)
    
# The results of the previous code cell indicate the largest component 
# contains at least 97% of the metabolites in the cell. I believe we can 
# safely discard the remainder. These nodes are also discarded from the digraph.

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
        print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
    
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


#%% Actual program operation
# Will go at end of file. Temporarily moved above working code to ensure all
# prereq functions have been called.    
tribeSampleDict =  importTaxonomy('../ExternalData/taxonomySAGs.csv')
createTribalGraph(tribeSampleDict)
computeGraphStats('../'+processedDataDir)
reduceToLargeComponent('../'+processedDataDir)