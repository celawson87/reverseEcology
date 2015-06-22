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

#%%

# Preliminaries
# This cell imports the packags and modules used throughout this notebook and
# creates the list of genomes to process. The user is allowed to define their
# local folder structure for data input and processing.

# Import Python packages.
# The packages 'cobra', 'networkx', 'numpy', and 'matplotlib' may not be 
# included with your Python distribution. They can be installed using pip. For
# more information, visit:
#   cobra: http://github.com/opencobra/cobrapy, used to process output from KBase
#   networkx: http://networkx.github.io/, used to perform graph analyses form which
# seed sets and scopes are computed
#   numpy: http://www.numpy.org/, adds support for matrices and arrays
#   matplotlib: http://matplotlib.org/, used for plotting
import cobra
import csv
import itertools
import os, glob
import networkx as nx

# Import scipy stack
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
import scipy.cluster.hierarchy as sch

# Import Python modules 
# These custom-written modules should have been included with the package
# distribution. They include functions for working with the SBML models from 
# KBase (sbmlFunctions and their graph representations (graphFunctions).
import sbmlFunctions as sf
import graphFunctions as gf

# Define data directories
# The user should ensure these directories match their data organization scheme
# processedDataDir: directory containing a folder for each genome. The folder
# name should match the genome name, and should contain a genome-scale model
# of the genome. The model should be named 'genomeNameBalanced.xml'.
# summaryStatsDir: directory for storing summary statistics
processedDataDir = 'ProcessedModelFiles'
summaryStatsDir = 'DataSummaries'
externalDataDir = 'ExternalData'

# Retrieve list of genomes to process by examing the contents of 
# 'processedDataDir', ignoring hidden folders. Subsequent computations are 
# performed by iterating over this list.
dirList =[]
for item in os.listdir('../'+processedDataDir):
    if not item.startswith('.'):
        dirList.append(item)

numSubDir = len(dirList)


#%%

# SBML to Adjacency List
# This cell converts each genome scale model from an SBML file to an adjacency
# list. An adjacency list is a collection of lists, one for each vertex in the
# graph, representing that vertex's connected edges. For more details, visit
# http://en.wikipedia.org/wiki/Adjacency_list
# Adjacency lists for each genome-scale model are written as text files in 
# each genome directory. Summary statistics about each graph are written in the
# summaryStatsDir as well.

# Create an array to store summary statistics. The array has three integer
# columns, which will contain the number of genes, metabolites, and reactions
# in the SBML file. 
modelStatArray = np.empty([numSubDir, 3], dtype = int)

# Create a file to record the summary statistics.
modelFile = open('../'+summaryStatsDir+'/'+'ModelStatistics.txt', 'w')
modelFile.write('Model,Genes,Metabolites,Reactions\n')

# Iterate over the list of genome directories. For each genome, read in the
# SBML file and update the 'description' field with the genome name. The number
# of genes, metabolites, and reactions in the SBML file is recorded in the
# 'modelStatArray' and written to 'modelFile.' Finally, the genome-scale model
# is converted to an adjacency list and written to file.
count = 0
print 'Converting SBML file to Adjacency List'

# Create an empty dictionary to store metabolite IDs and names
namesDict = {}

for curDir in dirList:
# Read in SBML file    
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
    model = cobra.io.read_sbml_model('../'+processedDataDir+'/'+curDir+'/'+curDir+'Balanced.xml')

# Create dictionary of metabolite names
    for metab in model.metabolites:
        namesDict[metab.id] = metab.name

# Update description field
    model.description = curDir;

# Read model statistics by invoking sbmlFunctions.getModelStats
    modelStatArray[count:] = sf.getModelStats(model)
    modelFile.write('%s,%i,%i,%i\n' % ('../'+processedDataDir+'/'+curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )
# Create adjacency list and write to file
    sf.adjacencyListFromModel(model, processedDataDir)
    count = count + 1

# Close files containing summary data
modelFile.close()

# Write completed dictionary to file as a csv file ExternalData/metabMap.csv
writer = csv.writer(open('../'+externalDataDir+'/'+'metabMap.csv', 'wb'))
for key, value in namesDict.items():
   writer.writerow([key, value])

#%%

# Statistics on Adjacency Lists
# This cell reads in the adjacency lists from the previous step and creates
# graph and directed graph (digraph) representations of each list. The objects
# are created using the networkX package. Summary statistics for the graph
# and directed graph are also reported and written to file.

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

#%%

# Reduction to largest component.

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


#%%

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
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir

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

# Plot summary statistics. The function plotSeedStats plots histograms of:
#   number of seed sets
#   size of seed sets
# The function plotSeedStats also plots scatter grams of the number of seed 
# sets vs the size of the original genome-scale model, either reactions or
# metabolites.
gf.plotSeedStats(seedSetList, reducedGraphStatArray, modelStatArray)