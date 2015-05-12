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

def getModelStats(model):
    statRow = [0]*3
    statRow[0] = len(model.genes)
    statRow[1] = len(model.metabolites)
    statRow[2] = len(model.reactions)
    return statRow
    
def adjacencyListFromModel(model):
# Establish a file to write the adjancey list to
    myFile = open(model.description+'/'+model.description+'AdjList.txt', 'w')
# For each reaction, loop over the reactants. For each reactant, loop over the 
# reaction products and construct the adjacency list for the reactant. If a 
# reaction is reversible, repeat the process in reverse.
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