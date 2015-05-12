###############################################################################
# sbmlToAdjacency.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Master function for performing reverse ecology simulations.
################################################################################

#%%
# Import python packages
import cobra
import itertools
import os, glob
import networkx as nx
import numpy as np
from matplotlib import pyplot

# Import additional functions which are found in relevant scripts
import sbmlFunctions as sf
import graphFunctions as gf

# Define data directories
processedDataDir = 'ProcessedModelFiles'
summaryStatsDir = 'DataSummaries'

dirList =[]
# Retrieve listing of model subdirectories
for item in os.listdir('../'+processedDataDir):
    if not item.startswith('.'):
        dirList.append(item)

numSubDir = len(dirList)


#%%
# Convert SBML model files to adjacency lists. Retreive model statistics.
modelStatArray = np.empty([numSubDir, 3], dtype = int)

modelFile = open('../'+summaryStatsDir+'/'+'ModelStatistics.txt', 'w')
modelFile.write('Model,Genes,Metabolites,Reactions\n')

count = 0
print 'Converting SBML file to Adjacency List'
for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
    model = cobra.io.read_sbml_model('../'+processedDataDir+'/'+curDir+'/'+curDir+'Balanced.xml')
    model.description = curDir;
# Read model statistics
    modelStatArray[count:] = sf.getModelStats(model)
    modelFile.write('%s,%i,%i,%i\n' % ('../'+processedDataDir+'/'+curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )
# Create adjacency list and write to file
    sf.adjacencyListFromModel(model, processedDataDir)
    count = count + 1
modelFile.close()


#%%
# Import adjacency lists as graphs. Retrieve graph statistics.
graphStatArray = np.empty([numSubDir, 4], dtype = int)
diGraphStatArray = np.empty([numSubDir, 4], dtype = int)

graphFile = open('../'+summaryStatsDir+'/'+'GraphStatistics.txt', 'w')
graphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

diGraphFile = open('../'+summaryStatsDir+'/'+'DiGraphStatistics.txt', 'w')
diGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

graphList = []
diGraphList = []

count = 0
print 'Computing Graph Statistics'
for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
# Read in adjacency list and convert to graph object
    myGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                              delimiter='\t', create_using=nx.Graph())
    myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())                            
    graphList.append(myGraph)
    diGraphList.append(myDiGraph)
# Read graph statistics                       
    graphStatArray[count:] = gf.getGraphStats(myGraph)
    graphFile.write('%s,%i,%i,%i,%i\n' % (curDir, graphStatArray[count,0], 
                                       graphStatArray[count,1], 
                                       graphStatArray[count, 2],
                                       graphStatArray[count, 3] ) )
    diGraphStatArray[count:] = gf.getDiGraphStats(myDiGraph)
    diGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, diGraphStatArray[count,0], 
                                       diGraphStatArray[count,1], 
                                       diGraphStatArray[count, 2],
                                       diGraphStatArray[count, 3] ) )

    count = count + 1
graphFile.close()
diGraphFile.close()

gf.plotGraphStats(graphStatArray)

#%%
# Reduce each graph to its largest component and write to file
# From the undirected graph, identify the nodes belonging to each component.
# Aggregate the remaining nodes and remove them from both the graph and the
# digraph.
# Write the diGraph to file as an adjacency list.

reducedGraphStatArray = np.empty([numSubDir, 4], dtype = int)
reducedDiGraphStatArray = np.empty([numSubDir, 4], dtype = int)

reducedGraphFile = open('../'+summaryStatsDir+'/'+'ReducedGraphStatistics.txt', 'w')
reducedGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

reducedDiGraphFile = open('../'+summaryStatsDir+'/'+'ReducedDiGraphStatistics.txt', 'w')
reducedDiGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

reducedGraphList = []
reducedDiGraphList = []
seedSetList = []

print 'Reducing to Largest Component'
count = 0

for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
    myGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                              delimiter='\t', create_using=nx.Graph())
    myDiGraph = nx.read_adjlist('../'+processedDataDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())                            
    
    subGraphs = sorted(nx.connected_components(myGraph), key = len, reverse=True)
    removeNodes = list(itertools.chain(*subGraphs[1:len(subGraphs)]))
#    print removeNodes
    
    myGraph.remove_nodes_from(removeNodes)
    myDiGraph.remove_nodes_from(removeNodes)
   
    reducedGraphList.append(myGraph)
    reducedDiGraphList.append(myDiGraph)

# Read graph statistics                       
    reducedGraphStatArray[count:] = gf.getGraphStats(myGraph)
    reducedGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, reducedGraphStatArray[count,0], 
                                       reducedGraphStatArray[count,1], 
                                       reducedGraphStatArray[count, 2],
                                       reducedGraphStatArray[count, 3] ) )
    reducedDiGraphStatArray[count:] = gf.getDiGraphStats(myDiGraph)
    reducedDiGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, reducedDiGraphStatArray[count,0], 
                                       reducedDiGraphStatArray[count,1], 
                                       reducedDiGraphStatArray[count, 2],
                                       reducedDiGraphStatArray[count, 3] ) )

# Write the graph as an adjancecy list
    nx.write_adjlist(myDiGraph, '../'+processedDataDir+'/'+curDir+'/'+curDir+'RedAdjList.txt')

# Compute seed sets and write to file
    mySCCList = list(nx.strongly_connected_components_recursive(myDiGraph))
    myCondensation = nx.condensation(myDiGraph)
    mySeeds = []    
    for node in myCondensation.nodes():
        inDeg = myCondensation.in_degree(node)
        if inDeg == 0:
            mySeeds.append(mySCCList[node])
    
    seedSetList.append(mySeeds)
    
    seedSets = open('../'+processedDataDir+'/'+curDir+'/'+curDir+'SeedSets.txt', 'w')
    for item in mySeeds:
        seedSets.write("%s\n" % item)
    seedSets.close()

    count = count + 1

reducedGraphFile.close()
reducedDiGraphFile.close()

gf.plotSeedStats(seedSetList, reducedGraphStatArray, modelStatArray)