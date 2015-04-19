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
import networkx as nx
import numpy as np
import os
import pandas
from matplotlib import pyplot

# Import additional functions which are found in relevant scripts
import sbmlFunctions as sf
import graphFunctions as gf

# Retrieve listing of model subdirectories
dirList = filter(os.path.isdir, glob.glob('*'))
numSubDir = len(dirList)

#%%
# Convert SBML model files to adjacencyLists. Retreive model statistics.
modelStatArray = np.empty([numSubDir, 3], dtype = int)

count = 0
for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
    model = cobra.io.read_sbml_model(curDir+'/'+curDir+'Balanced.xml')
    model.description = curDir;
# Read model statistics
    modelStatArray[count:] = sf.getModelStats(model)
    sbmlFunctions.adjacencyListFromModel(model)
    count = count + 1

# Write model statistics to file
colLabels = ['Genes', 'Metabolites', 'Reactions']
modelStatDF = pandas.DataFrame(modelStatArray, columns=colLabels, index=dirList) 
modelStatDF.to_csv('ModelStatistics.txt')

#%%
# Import adjacencylists as graphs. Retrieve graph statistics.
graphStatArray = np.empty([numSubDir, 4], dtype = int)

count = 0
for curDir in dirList:
    print 'Processing directory', count+1, 'of', numSubDir, ':', curDir
# Construct adjacency list and convert to graph object
    myGraph = nx.read_adjlist(curDir+'/'+curDir+'AdjList.txt', delimiter='\t', 
                              create_using=nx.Graph())
    myDiGraph = nx.read_adjlist(curDir+'/'+curDir+'AdjList.txt', delimiter='\t', 
                              create_using=nx.DiGraph())
# Read graph statistics                       
    graphStatArray[count:] = gf.getGraphStats(myGraph)
    count = count + 1

colLabels = ['Nodes', 'Edges', 'TotalComponents', 'Size of Largest']
graphStatDF = pandas.DataFrame(graphStatArray, columns=colLabels, index=dirList) 
graphStatDF.to_csv('GraphStatistics.txt')

# Make some graphs
gf.plotGraphStats(graphStatArray)

#%%

# To identify the most highly connected metabolites we will:
curDir = 'TBhypometabat3838' #smallest
#curDir - 'TBhypometabat3815' #largest
myGraph = nx.read_adjlist(curDir+'/'+curDir+'AdjList.txt', delimiter='\t', 
                              create_using=nx.Graph())
myCent = nx.degree_centrality(myGraph)

import operator
myCentSorted = sorted(myCent.items(), key=operator.itemgetter(1), reverse=True)
myRank = np.linspace(1,len(myCentSorted), len(myCentSorted))

# Plot the data
pyplot.scatter(myRank, myCentSorted)
# Not plotting properly b/c myCentSorted is a list of lists. How do I extract
# the second column to plot it?


# Loop over all graphs. For each node, check if node has already been found. If
# it has a higher connectivity, then update the hash.
