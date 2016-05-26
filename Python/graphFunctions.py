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
import pandas as pd
import re
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
    binWidth = 10
    binMin = (min(graphStatArray[:,0]) / binWidth)*binWidth
    binMax = ((max(graphStatArray[:,0]) / binWidth)*binWidth) + binWidth
    plt.hist(graphStatArray[:,0], bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
    plt.xlabel('Total Nodes')
    plt.xlim(binMin, binMax)
    plt.ylabel('Fraction of Graphs')

# Histogram of number of components
    myWeight = np.ones_like(graphStatArray[:,2]) / float(len(graphStatArray[:,2]))
    plt.figure(2)
    binWidth = 1
    binMin = min(graphStatArray[:,2])
    binMax = max(graphStatArray[:,2]) + binWidth
    plt.hist(graphStatArray[:,2], bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
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
    binWidth = 10
    binMin = (min(myNumSeedSets) / binWidth)*binWidth
    binMax = ((max(myNumSeedSets) / binWidth)*binWidth) + binWidth
    plt.hist(myNumSeedSets, bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
    plt.xlabel('Number of Seed Sets')
    plt.ylabel('Fraction of Graphs')
 
# Histogram of size of individual seed sets
    myWeight = np.ones_like(mySizeOfSeedSets) / float(len(mySizeOfSeedSets))
    plt.figure(8)
    binWidth = 1
    binMin = min(mySizeOfSeedSets)
    binMax = max(mySizeOfSeedSets) + binWidth
    plt.hist(mySizeOfSeedSets, bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
    plt.xlabel('Metabolites in Seed Set')
    plt.xlim(0, max(mySizeOfSeedSets))
    plt.ylim(0, 1)
    plt.ylabel('Fraction of Seed Sets')
    
# Zoomed histogram of size of individual seed sets
    plt.figure(9)
    [n, bins, patches] = plt.hist(mySizeOfSeedSets, bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
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

def plotSeedStatsForTribes(seedSetList, reducedGraphStatArray):
# Histogram of total number of seed sets
    myNumSeedSets = []
    mySizeOfSeedSets = []
    
    for setOfSeedSets in seedSetList:
        myNumSeedSets.append(len(setOfSeedSets))
        for seedSet in setOfSeedSets:
            mySizeOfSeedSets.append(len(seedSet))

    myWeight = np.ones_like(myNumSeedSets) / float(len(myNumSeedSets))
    plt.figure(12)
    binWidth = 10
    binMin = (min(myNumSeedSets) / binWidth)*binWidth
    binMax = ((max(myNumSeedSets) / binWidth)*binWidth) + binWidth
    plt.hist(myNumSeedSets, bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
    plt.xlabel('Number of Seed Sets')
    plt.ylabel('Fraction of Graphs')
 
# Histogram of size of individual seed sets
    myWeight = np.ones_like(mySizeOfSeedSets) / float(len(mySizeOfSeedSets))
    plt.figure(13)
    binWidth = 1
    binMin = min(mySizeOfSeedSets)
    binMax = max(mySizeOfSeedSets) + binWidth
    plt.hist(mySizeOfSeedSets, bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
    plt.xlabel('Metabolites in Seed Set')
    plt.xlim(0, max(mySizeOfSeedSets)+1)
    plt.ylim(0, 1)
    plt.ylabel('Fraction of Seed Sets')
    
# Zoomed histogram of size of individual seed sets
    plt.figure(14)
    [n, bins, patches] = plt.hist(mySizeOfSeedSets, bins=range(binMin, binMax + binWidth, binWidth), weights=myWeight)
    plt.xlabel('Metabolites in Seed Set (Zoomed)')
    plt.xlim(0, max(mySizeOfSeedSets)+1)
    plt.ylim(0, 1.1*n[1])
    plt.ylabel('Fraction of Seed Sets (Zoomed)')
    
    return
    
################################################################################
    
# createMergedGraph
# In this function, all samples from a tribe are identified. Each sample is
# converted to a graph object and merged with the previous graph. The final
# graph is written to file.

def createMergedGraph(groupSampleDict, processedDataDir, rawModelDir):

    print 'Merging genomes from specified taxonomic groups (lineage/clade/group)'
    
# Loop over the keys of the dictionary, one for each group
    for group in groupSampleDict:

# Create an empty graph object
        mergedGraph = nx.DiGraph()

# Read in the graph of the group and merge with the graph from the previous
# iteration
        for sample in groupSampleDict[group]:

# Read in adjacency list and convert to digraph object
            myDiGraph = nx.read_adjlist('../'+rawModelDir+'/'+sample+'/'+sample+'AdjList.txt',
                                create_using=nx.DiGraph())

# Append to the previous graph
            mergedGraph = nx.compose(mergedGraph, myDiGraph)

# Check that the proper output directory exists. It not, create it.
        if not os.path.exists('../'+processedDataDir+'/'+group):
            os.makedirs('../'+processedDataDir+'/'+group)
    
        nx.write_adjlist(mergedGraph, '../'+processedDataDir+'/'+group+'/'+group+'AdjList.txt')
        nx.write_graphml(mergedGraph, '../'+processedDataDir+'/'+group+'/'+group+'Graph.xml')

    return

################################################################################
    
# computeGraphStats
# This functions reads in the adjacency lists from the given directory and 
# creates graph and directed graph (digraph) representations of each list. The 
# objects are created using the networkX package. Summary statistics for the 
# graph and directed graph are also reported and written to file.

def computeGraphStats(dirList, processedDataDir, summaryStatsDir):
    
# Check that folders exist and create them if necessary
    if not os.path.exists('../'+summaryStatsDir):
        os.makedirs('../'+summaryStatsDir)
        
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
                              create_using=nx.Graph())

# Read in adjacency list and convert to digraph object
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                create_using=nx.DiGraph())                            

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
                              create_using=nx.Graph())

# Read in adjacency list and convert to digraph object
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                create_using=nx.DiGraph())                            

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
        nx.write_graphml(myDiGraph, '../'+processedDataDir+'/'+curDir+'/'+curDir+'RedGraph.xml')
                                       
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
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'RedAdjList.txt',
                                create_using=nx.DiGraph())                            
    
    # Compute the list of SCCs for the digraph as well as its condensation
        myCondensation = nx.condensation(myDiGraph)
        nx.write_adjlist(myCondensation, '../'+processedDataDir+'/'+curDir+'/'+curDir+'SCCAdjList.txt')
    
    # For some reason, the condensation cannot be written to GraphML. Instead, re-read the 
    # adjacency list and write that to GraphML.
        myTempGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'SCCAdjList.txt',
                                create_using=nx.DiGraph())                            
        nx.write_graphml(myTempGraph, '../'+processedDataDir+'/'+curDir+'/'+curDir+'SCCGraph.xml')
    
    
    # Invert the mapping dictionary to map SCC nodes to their original compoundsm
        mapDict = dict()
        for key in myCondensation.graph.items()[0][1].keys():
            value = str(myCondensation.graph.items()[0][1][key])
        # If the value exists as a key in mapDict, append the new value
            if value in mapDict.keys():
                mapDict[value].append(str(key))
        # Otherwise create it
            else:
                mapDict[value] = [str(key)]
    
        dictFile=open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SCCDict.txt', "w")
        for key in mapDict.keys():
            dictFile.write(str(key)+',')
            dictFile.write(",".join(str(value) for value in mapDict[key]))
            dictFile.write('\n')
        dictFile.close()
        
    # "List of lists" of seed metabolites. Each element is a list of nodes belonging
    # to an SCC which is also a seed set.
        mySeeds = []    
    
    # For each node (SCC) of the condensation, examine each its in-degree. If the
    # in-degree is zero (only outgoing edges), the SCC is a seed set. Append the
    # SCC (as a list of nodes) to the list of seed sets.
        for node in myCondensation.nodes():
            inDeg = myCondensation.in_degree(node)
            if inDeg == 0:
                mySeeds.append(mapDict[str(node)])
        seedSetList.append(mySeeds)
    
    # Update the list of seed metabolites: replace the Model SEED metabolite 
    # identifier with its common name. Note: The file metabMap.csv was created 
    # manually from the seed database, and should be updated to reflect the 
    # particulars of your data set. As above, record the seed metabolite for each
    # graph.
    
    # First read metabMap.csv in as a dictionary
        with open('../'+externalDataDir+'/'+'metabMap.csv', mode='rU') as inFile:
            reader = csv.reader(inFile)
            namesDict = dict((rows[0],rows[1]) for rows in reader)
    
    # Compute weights for each seed metabolite and write to file. Each row of the 
    # output file contains a metabolite and its weight (1 / size of the seed set). 
    
        seedFile = open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedCompounds.txt', 'w')
        for seed in mySeeds:
            myWeight = 1 / float(len(seed))
            for metab in seed:
                seedFile.write('%s\t%s\t%f\n' % (metab, namesDict[re.sub('_[a-d]', '', metab)], myWeight) )
        seedFile.close()
    
        count = count + 1
    
    return seedSetList
    
################################################################################
# definition of Metabolite Groups
# defineMetabGroups

# The acI metabolic networks have a bow-tie structure, with a single giant
# connected component (GCC) containing a majority of the metabolites. We 
# decompose the metabolites in the network into groups, depending on their 
# relationship to the GCC.

# Seed compounds which point to the GCC
# The remaining seed compounds
# Sink compounds (no outward arcs) which come from the GCC
# The remaining sink compounds
# All other compounds

def defineMetabGroups(processedDataDir, level, taxonFile):

    # Obtain the groupList
    taxonClass = pd.DataFrame.from_csv(taxonFile, sep=',')
    taxonClass = taxonClass.dropna()
    
    # Extract the unique tribes found in the dataset
    if level=='Genome':
        groupList = list(taxonClass.index.values)
    else:   
        groupList = pd.unique(taxonClass[level].values)
        groupList.sort(axis=0)
        groupList = [ group for group in groupList if not group.startswith('Unknown') ]
        groupList = sorted(groupList, key=str.lower)
        
    for curDir in groupList:
        # Read in the clade-level network and compute its condensation
        myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'RedAdjList.txt',
                                        create_using=nx.DiGraph())                            
    
        # Compute the list of SCCs for the digraph as well as its condensation
        myCondensation = nx.condensation(myDiGraph)    
       
       # Invert the mapping dictionary to map SCC nodes to their original compoundsm
        mapDict = dict()
        for key in myCondensation.graph.items()[0][1].keys():
            value = str(myCondensation.graph.items()[0][1][key])
            # If the value exists as a key in mapDict, append the new value
            if value in mapDict.keys():
                mapDict[value].append(str(key))
                # Otherwise create it
            else:
                mapDict[value] = [str(key)]
    
        # Compute the seed compounds                    
        seedSetList = []
        sccSeedSetList = []
        GCC = ''
        for node in myCondensation.nodes():
            inDeg = myCondensation.in_degree(node)
            # Check the in-degree. If 0, it's a seed.
            if inDeg == 0:
                seedSetList.append(mapDict[str(node)])
                sccSeedSetList.append(node)
            # The GCC has the highest in-degree. If the node has a higher in-degree,
            # is is the new GCC
            if len(GCC) == 0:
                GCC = str(node)
            if inDeg > myCondensation.in_degree(int(GCC)):
                GCC = str(node)
        # Flatten the list
        seedSetList = [seed for seedList in seedSetList for seed in seedList]
        
        # Identify the seed compounds which point to the GCC
        # Work with seeds in the SCC graph, then convert to compounds
        sccToGccList = []
        seedToGccList = []
        
        for node in sccSeedSetList:
            if int(GCC) in nx.all_neighbors(myCondensation, node):
                sccToGccList.append(node)
        for node in sccToGccList:
            seedToGccList.append(mapDict[str(node)])
        seedToGccList = [seed for seedList in seedToGccList for seed in seedList]
        
        # And those that don't
        remainingSeedList = [seed for seed in seedSetList if seed not in seedToGccList]
        
        # Write these lists to file
        with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'seedsToGCC.txt', "w") as outFile:
            for seed in seedToGccList:
                outFile.write(seed+'\n')
        with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'seedsToNotGCC.txt', "w") as outFile:
            for seed in remainingSeedList:
                outFile.write(seed+'\n')
    
        # Compute the sink compounds                    
        sinkSetList = []
        sccSinkSetList = []
        for node in myCondensation.nodes():
            outDeg = myCondensation.out_degree(node)
            # Check the in-degree. If 0, it's a seed.
            if outDeg == 0:
                sinkSetList.append(mapDict[str(node)])
                sccSinkSetList.append(node)
            # The GCC has the highest in-degree. If the node has a higher in-degree,
            # is is the new GCC
        # Flatten the list
        sinkSetList = [sink for sinkList in sinkSetList for sink in sinkList]
        
        # Identify the seed compounds which point to the GCC
        # Work with seeds in the SCC graph, then convert to compounds
        sccToGccList = []
        sinkToGccList = []
        
        for node in sccSinkSetList:
            if int(GCC) in nx.all_neighbors(myCondensation, node):
                sccToGccList.append(node)
        for node in sccToGccList:
            sinkToGccList.append(mapDict[str(node)])
        sinkToGccList = [sink for sinkList in sinkToGccList for sink in sinkList]
        
        # And those that don't
        remainingSinkList = [sink for sink in sinkSetList if sink not in sinkToGccList]
        
        # Write these lists to file
        with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'sinkFromGCC.txt', "w") as outFile:
            for sink in sinkToGccList:
                outFile.write(sink+'\n')
        with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'sinkFromNotGCC.txt', "w") as outFile:
            for sink in remainingSinkList:
                outFile.write(sink+'\n')
    
        # Identify the remaining metabolites and write to file
        remainingMetabList = [str(metab) for metab in myDiGraph.nodes()]
        remainingMetabList = [metab for metab in remainingMetabList if metab not in seedToGccList]
        remainingMetabList = [metab for metab in remainingMetabList if metab not in remainingSeedList]
        remainingMetabList = [metab for metab in remainingMetabList if metab not in sinkToGccList]
        remainingMetabList = [metab for metab in remainingMetabList if metab not in remainingSinkList]
        
        with open('../'+processedDataDir+'/'+curDir+'/'+curDir+'remainingMetabs.txt', "w") as outFile:
            for metab in remainingMetabList:
                outFile.write(metab+'\n')
    return