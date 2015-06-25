###############################################################################
# graphFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with networkx graph objects.
################################################################################

# Import Python packages.
#import math
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
#from collections import Counter

import csv
import os
import itertools

################################################################################

# getGraphStats
# Function to retrieve statistics about a graph
# Input: network object of a graph
# Output: array with four integer columns, containing the number of nodes 
# (metabolites), edges, total components, and size of the largest component.

def getGraphStats(graph):
    statRow = [0]*4
    statRow[0] = graph.number_of_nodes()
    statRow[1] = graph.number_of_edges()
    myConComp = nx.connected_components(graph)
    myConCompList = sorted(myConComp, key = len, reverse=True)
    statRow[2] = len(myConCompList)
    statRow[3] = len(myConCompList[0])
    return statRow

################################################################################

# getDiGraphStats
# Function to retrieve statistics about a digraph
# Input: network object of a digraph
# Output: array with four integer columns, containing the number of nodes 
# (metabolites), edges, total components, and size of the largest component.

def getDiGraphStats(diGraph):
    statRow = [0]*4
    statRow[0] = diGraph.number_of_nodes()
    statRow[1] = diGraph.number_of_edges()
    myConComp = nx.strongly_connected_components(diGraph)
    myConCompList = sorted(myConComp, key = len, reverse=True)
    statRow[2] = len(myConCompList)
    statRow[3] = len(myConCompList[0])
    return statRow

################################################################################

# plotGraphStats
# Plot summary statistics of a collection of graph objects. The function plots
# historams of:
#   graph size (number of nodes)
#   total number of components
#   size of largest compmonent, as fraction of total nodes
# Input: array containing one row for each graph object. Array columns 
# correspond to: the number of nodes (metabolites), edges, total components, 
# and size of the largest component.
# Output: collection of plots

def plotGraphStats(graphStatArray):
# Histogram of number of nodes
    myWeight = np.ones_like(graphStatArray[:,0]) / float(len(graphStatArray[:,0]))
    plt.figure(1)
    plt.hist(graphStatArray[:,0], weights=myWeight)
    plt.xlabel('Total Nodes')
    plt.ylabel('Fraction of Graphs')

# Histogram of number of components
    myWeight = np.ones_like(graphStatArray[:,2]) / float(len(graphStatArray[:,2]))
    plt.figure(2)
    plt.hist(graphStatArray[:,2], weights=myWeight)
    plt.xlabel('Number of Components')
    plt.ylabel('Fraction of Graphs')

# Histogram of largest component size
    myArray = np.true_divide(graphStatArray[:,3], graphStatArray[:,0]);
    myWeight = np.ones_like(myArray) / float(len(myArray))
    plt.figure(3)
    plt.hist(myArray, weights=myWeight)
    plt.xlabel('Fraction of Nodes in Largest Component')
    plt.ylabel('Fraction of Graphs')
    
    return

################################################################################

# plotDiGraphStats
# Plot summary statistics of a collection of digraph objects. The function plots
# historams of:
#   digraph size (number of nodes)
#   total number of components
#   size of largest compmonent, as fraction of total nodes
# Input: array containing one row for each digraph object. Array columns 
# correspond to: the number of nodes (metabolites), edges, total components, 
# and size of the largest component.
# Output: collection of plots

def plotDiGraphStats(diGraphStatArray):
# Histogram of number of nodes
    myWeight = np.ones_like(diGraphStatArray[:,0]) / float(len(diGraphStatArray[:,0]))
    plt.figure(4)
    plt.hist(diGraphStatArray[:,0], weights=myWeight)
    plt.xlabel('Total Nodes')
    plt.ylabel('Fraction of Graphs')

# Histogram of number of components
    myWeight = np.ones_like(diGraphStatArray[:,2]) / float(len(diGraphStatArray[:,2]))
    plt.figure(5)
    plt.hist(diGraphStatArray[:,2], weights=myWeight)
    plt.xlabel('Number of Seed Sets')
    plt.ylabel('Fraction of Genomes')

# Histogram of largest component size
    myWeight = np.ones_like(diGraphStatArray[:,3]) / float(len(diGraphStatArray[:,3]))
    plt.figure(6)
    plt.hist(diGraphStatArray[:,3], weights=myWeight)
    plt.xlabel('Number of Compounds in Largest Seed Set')
    plt.ylabel('Fraction of Genomes')

    return

################################################################################

# plotSeedStats
# Plot summary statistics of a collection of seed sets. The function plots
# historams of:
#   digraph size (number of nodes)
#   total number of components
#   size of largest compmonent, as fraction of total nodes
# Inputs:
#  seedSetList: "List of lists" of seed metabolites. Each element is a list of nodes belonging
#   to an SCC which is also a seed set.
#  reducedGraphStatArray: array containing one row for each graph object. Array 
#   columns correspond to: the number of nodes (metabolites), edges, total 
#   components, and size of the largest component.
#  modelStatArray: array containing one row for each original SBML file/model. 
#    Array columns correspond to:  the number of genes, metabolites, and reactions
#    in the SBML file. 
# Output: collection of plots

def plotSeedStats(seedSetList, reducedGraphStatArray, modelStatArray):
# Histogram of total number of seed sets
    myNumSeedSets = []
    mySizeOfSeedSets = []
    
    for setOfSeedSets in seedSetList:
        myNumSeedSets.append(len(setOfSeedSets))
        for seedSet in setOfSeedSets:
            mySizeOfSeedSets.append(len(seedSet))
        
    myWeight = np.ones_like(myNumSeedSets) / float(len(myNumSeedSets))
    plt.figure(7)
    plt.hist(myNumSeedSets, weights=myWeight)
    plt.xlabel('Number of Seed Sets')
    plt.ylabel('Fraction of Graphs')
 
# Histogram of size of individual seed sets
    myWeight = np.ones_like(mySizeOfSeedSets) / float(len(mySizeOfSeedSets))
    plt.figure(8)
    plt.hist(mySizeOfSeedSets, weights=myWeight)
    plt.xlabel('Metabolites in Seed Set')
    plt.xlim(0, max(mySizeOfSeedSets))
    plt.ylim(0, 1)
    plt.ylabel('Fraction of Seed Sets')
    
# Zoomed histogram of size of individual seed sets
    plt.figure(9)
    [n, bins, patches] = plt.hist(mySizeOfSeedSets, weights=myWeight)
    plt.xlabel('Metabolites in Seed Set (Zoomed)')
    plt.xlim(0, max(mySizeOfSeedSets))
    plt.ylim(0, 1.1*n[1])
    plt.ylabel('Fraction of Seed Sets (Zoomed)')
    
# Scatter plots of seed sets ize compared to model size, as measured by
# metabolites or reactions
    plt.figure(10)
    plt.scatter(reducedGraphStatArray[:,0], myNumSeedSets)
    plt.xlabel('Total Number of Metabolites')
    plt.ylabel('Number of Seed Sets')
    
    plt.figure(11)
    plt.scatter(modelStatArray[:,2], myNumSeedSets)
    plt.xlabel('Total Number of Reactions')
    plt.ylabel('Number of Seed Sets')
    
    return
    
################################################################################
    
# createTribalGraph
# In this function, all samples from a tribe are identified. Each sample is
# converted to a graph object and merged with the previous graph. The final
# graph is written to file.

def createTribalGraph(tribeSampleDict, processedDataDir, rawModelDir):

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

# Check that the proper output directory exists. It not, create it.
        if not os.path.exists('../'+processedDataDir+'/'+tribe):
            os.makedirs('../'+processedDataDir+'/'+tribe)
    
        nx.write_adjlist(tribalGraph, '../'+processedDataDir+'/'+tribe+'/'+tribe+'AdjList.txt', delimiter='\t')

    return

################################################################################
    
# computeGraphStats
# This functions reads in the adjacency lists from the given directory and 
# creates graph and directed graph (digraph) representations of each list. The 
# objects are created using the networkX package. Summary statistics for the 
# graph and directed graph are also reported and written to file.

def computeGraphStats(dirList, processedDataDir, summaryStatsDir):
    
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
    print 'Computing graph statistics'

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
        graphStatArray[count:] = getGraphStats(myGraph)
        graphFile.write('%s,%i,%i,%i,%i\n' % (curDir, graphStatArray[count,0], 
                                       graphStatArray[count,1], 
                                       graphStatArray[count, 2],
                                       graphStatArray[count, 3] ) )

# Read model statistics by invoking graphFunctions.getDiGraphStats
        diGraphStatArray[count:] = getDiGraphStats(myDiGraph)
        diGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, diGraphStatArray[count,0], 
                                       diGraphStatArray[count,1], 
                                       diGraphStatArray[count, 2],
                                       diGraphStatArray[count, 3] ) )

        count = count + 1
    
# Close files containing summary data
    graphFile.close()
    diGraphFile.close()
    
    return graphStatArray, diGraphStatArray
    
################################################################################

# reduceToLargeComponent
# This function iterates over a list of genomes and identifies the largest
# component of that genome's network graph. Nodes outside of this component are
# discarded, and the reduced graph is written to file.
    
def reduceToLargeComponent(dirList, processedDataDir, summaryStatsDir):
    
    numSubDir = len(dirList)

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
    print 'Reducing to largest component'

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
        reducedGraphStatArray[count:] = getGraphStats(myGraph)
        reducedGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, reducedGraphStatArray[count,0], 
                                       reducedGraphStatArray[count,1], 
                                       reducedGraphStatArray[count, 2],
                                       reducedGraphStatArray[count, 3] ) )
# Read model statistics by invoking graphFunctions.getDiGraphStats
        reducedDiGraphStatArray[count:] = getDiGraphStats(myDiGraph)
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
    
    return reducedGraphStatArray
    
################################################################################

# Computation of seed sets.
# computeSeedSets

# This function computes the seed compounds of a metabolic network graph.

# First, it computes the strongly connected components (SCCs) of the reduced 
# digraph (from function 'reduceToLargestComponent'). An SCC is a group of
# nodes, such that from each node there exists a path to all other nodes in 
# the component. SCCs are candidates for seed sets.

# Second, SCCs are evaluated to see if they are seed sets: any SCC with no 
# outgoing edges is a seed set. This is done by converting the digraph to its 
# condensation (a new graph) in which each SCC is represented as a single
# node. 

# Third, seed sets are written to file and summary statistics are computed
# for each seed set. Additional statistics on the reduced graph and digraph
# are also computed.
    
def computeSeedSets(dirList, externalDataDir, processedDataDir):
        
# Create lists to store seed sets
# seedSetList is a list of lists. Each outer list contains all the seed sets
# for that graph.
    seedSetList = []

# Iterate over the list of genome directories. For each reduced digraph, 
# identify its condensation (SCCs). For each node of the SCC, check if it
# is a seed set by computing its in-degree. If yes, append the SCC (as a list
# of nodes) to the list of seed sets. Then compute some summary statistics.
    count = 0
    print 'Computing seed sets'

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

    return seedSetList