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

# This module contains functions for working with networkx graph objects.

# Import Python packages.
# The packages 'networkx', 'numpy', and 'matplotlib' may not be 
# included with your Python distribution. They can be installed using pip. For
# more information, visit:
# networkx: http://networkx.github.io/, used to perform graph analyses form which
# seed sets and scopes are computed
# numpy: http://www.numpy.org/, adds support for matrices and arrays
# matplotlib: http://matplotlib.org/, used for plotting
import math
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

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