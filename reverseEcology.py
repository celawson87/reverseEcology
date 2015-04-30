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
import os, glob
import networkx as nx
import numpy as np
from matplotlib import pyplot

# Import additional functions which are found in relevant scripts
import sbmlFunctions as sf
import graphFunctions as gf

# Retrieve listing of model subdirectories
dirList = filter(os.path.isdir, glob.glob('*'))
numSubDir = len(dirList)


#%%
# Convert SBML model files to adjacency lists. Retreive model statistics.
modelStatArray = np.empty([numSubDir, 3], dtype = int)

count = 0
modelFile = open('ModelStatistics.txt', 'w')
modelFile.write('Model,Genes,Metabolites,Reactions\n')

for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
    model = cobra.io.read_sbml_model(curDir+'/'+curDir+'Balanced.xml')
    model.description = curDir;
# Read model statistics
    modelStatArray[count:] = sf.getModelStats(model)
    modelFile.write('%s,%i,%i,%i\n' % (curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )
# Create adjacency list and write to file
    sf.adjacencyListFromModel(model)
    count = count + 1
modelFile.close()


#%%
# Import adjacency lists as graphs. Retrieve graph statistics.

graphStatArray = np.empty([numSubDir, 4], dtype = int)

count = 0
graphFile = open('GraphStatistics.txt', 'w')
graphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

graphList = []
diGraphList = []
for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
# Read in adjacency list and convert to graph object
    myGraph = nx.read_adjlist(curDir+'/'+curDir+'AdjList.txt',
                              delimiter='\t', create_using=nx.Graph())
    myDiGraph = nx.read_adjlist(curDir+'/'+curDir+'AdjList.txt',
                                delimiter='\t', create_using=nx.DiGraph())                            
    graphList.append(myGraph)
    diGraphList.append(myDiGraph)
# Read graph statistics                       
    graphStatArray[count:] = gf.getGraphStats(myGraph)
    graphFile.write('%s,%i,%i,%i,%i\n' % (curDir, graphStatArray[count,0], 
                                       graphStatArray[count,1], 
                                       graphStatArray[count, 2],
                                       graphStatArray[count, 3] ) )
    count = count + 1
graphFile.close()