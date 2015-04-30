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

import networkx as nx
import numpy as np
from matplotlib import pyplot

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
    myWeight = np.ones_like(graphStatArray[:,0]) / 
                            float(len(graphStatArray[:,0]))
    pyplot.figure(1)
    pyplot.hist(graphStatArray[:,0], weights=myWeight)
    pyplot.xlabel('Total Nodes')
    pyplot.ylabel('Fraction of Graphs')

# Histogram of number of components
    myWeight = np.ones_like(graphStatArray[:,2]) / 
                            float(len(graphStatArray[:,2]))
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