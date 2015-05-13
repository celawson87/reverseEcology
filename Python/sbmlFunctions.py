###############################################################################
# sbmlFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for manipulating SBML files
################################################################################

# This module contains functions for working with cobrapy model objects built
# from SBML files.

# Function to retrieve statistics about a model
# Input: cobrapy model object of a genome-scale model
# Output: array with three integer columns, containing the number of genes, 
# metabolites, and reactions in the model
def getModelStats(model):
    statRow = [0]*3
    statRow[0] = len(model.genes)
    statRow[1] = len(model.metabolites)
    statRow[2] = len(model.reactions)
    return statRow

# Function to convert a cobrapy model object to an adjaceny list. Also writes
# adjacency list to file. For each reaction in the model, creates an edge
# between all (reactant, product) pairs. If a reaction is reversible, also
# creates edges between all (product, reactant) pairs.
# Input: cobrapy model object, model directory
# Output: None.
def adjacencyListFromModel(model, processedDataDir):
# Establish a file for the adjacency list
    myFile = open('../'+processedDataDir+'/'+model.description+'/'+model.description+'AdjList.txt', 'w')
# For each reaction, loop over the reactants. For each reactant, loop over the 
# reaction products and create an edge between the reactant and products. If a 
# reaction is reversible, repeat the process in reverse, creating an edge
# between each product and reactant.
    for myRxn in model.reactions:
        for myReactant in myRxn.reactants:
            myFile.write(myReactant.id+'\t')
            for myProduct in myRxn.products:
                myFile.write(myProduct.id+'\t')
            myFile.write('\n')
        if myRxn.reversibility == True:
            for myProduct in myRxn.products:
                myFile.write(myProduct.id+'\t')
                for myReactant in myRxn.reactants:
                    myFile.write(myReactant.id+'\t')
                myFile.write('\n')
    myFile.close()
    return