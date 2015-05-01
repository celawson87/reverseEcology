###############################################################################
# graphFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for manipulating SBML files
################################################################################

import math
import networkx as nx
import numpy as np
from matplotlib import pyplot
from collections import Counter


def getGraphStats(graph):
    statRow = [0]*4
    statRow[0] = graph.number_of_nodes()
    statRow[1] = graph.number_of_edges()
    myConComp = nx.connected_components(graph)
    myConCompList = sorted(myConComp, key = len, reverse=True)
    statRow[2] = len(myConCompList)
    statRow[3] = len(myConCompList[0])
    return statRow
    
def plotGraphStats(graphStatArray):
    myWeight = np.ones_like(graphStatArray[:,0]) / float(len(graphStatArray[:,0]))
    pyplot.figure(1)
    pyplot.hist(graphStatArray[:,0], weights=myWeight)
    pyplot.xlabel('Total Nodes')
    pyplot.ylabel('Fraction of Graphs')

# Histogram of number of components
    myWeight = np.ones_like(graphStatArray[:,2]) / float(len(graphStatArray[:,2]))
    pyplot.figure(2)
    pyplot.hist(graphStatArray[:,2], weights=myWeight)
    pyplot.xlabel('Number of Components')
    pyplot.ylabel('Fraction of Graphs')

# Histogram of largest component size
    myArray = np.true_divide(graphStatArray[:,3], graphStatArray[:,0]);
    myWeight = np.ones_like(myArray) / float(len(myArray))
    pyplot.figure(3)
    pyplot.hist(myArray, weights=myWeight)
    pyplot.xlabel('Fraction of Nodes in Largest Component')
    pyplot.ylabel('Fraction of Graphs')
    
    return
    
    
def findTopMetab(myPct, graphList):

    aggNodeCount = Counter()
    for graph in graphList:
        aggNodeCount = aggNodeCount + Counter(nx.degree(graph))
    
    aggNodeList = aggNodeCount.most_common()
    # Indices for splitting metabolites into two sets.
    totalNodes = len(aggNodeCount)
    pctIndex = int(math.ceil(totalNodes*myPct))
    
    # Plot the number of edges associated with each metabolite
    x = np.linspace(1, totalNodes, totalNodes)
    y = zip(*aggNodeList)[1]
    x0 = x[0:pctIndex-1]
    x1 = x[pctIndex:totalNodes]
    y0 = y[0:pctIndex-1]
    y1 = y[pctIndex:totalNodes]
    
    pyplot.loglog(x0, y0, marker='.', color='red', linestyle='none')
    pyplot.loglog(x1, y1, marker='.', color='black', linestyle='none')
    pyplot.xlim(0, len(x))
    pyplot.ylim(0, max(y))
    pyplot.xlabel('Metabolite Rank')
    pyplot.ylabel('Number of Edges')
    
    # Return a list of the most-connected metabolites
    metabFile = open('connectedMetabs.txt', 'w')
    metabFile.write('Metabolite, Edges\n')
    
    removeMetabs = []
    for count in range(pctIndex):
        metabFile.write('%s,%i\n' % (aggNodeList[count][0], 
                                     aggNodeList[count][1] ) )
        removeMetabs.append(aggNodeList[count][0])
                                 
    metabFile.close()

    return removeMetabs