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

# Import python modules
import cobra
import csv
import numpy as np

################################################################################
    
# dirListToAdjacencyList
# This function iterates over a list of genome directories and converts each
# genome scale model from an SBML file to an adjacency list. Adjacency lists 
# for each genome-scale model are written as text files in each genome 
# directory. Summary statistics about each graph are written in the
# summaryStatsDir as well.

def dirListToAdjacencyList(dirList, externalDataDir, processedDataDir, summaryStatsDir):

    numSubDir = len(dirList)

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
    print 'Converting SBML files to adjacency lists'

# Create an empty dictionary to store metabolite IDs and names
    namesDict = {}
    
    for curDir in dirList:
# Read in SBML file    
        model = cobra.io.read_sbml_model('../'+processedDataDir+'/'+curDir+'/'+curDir+'Balanced.xml')

# Create dictionary of metabolite names
        for metab in model.metabolites:
            namesDict[metab.id] = metab.name

# Update description field
        model.description = curDir;

# Read model statistics by invoking sbmlFunctions.getModelStats
        modelStatArray[count:] = getModelStats(model)
        modelFile.write('%s,%i,%i,%i\n' % ('../'+processedDataDir+'/'+curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )

# Create adjacency list and write to file
        adjacencyListFromModel(model, processedDataDir)
        count = count + 1

# Close files containing summary data
    modelFile.close()

# Write completed dictionary to file as a csv file ExternalData/metabMap.csv
    writer = csv.writer(open('../'+externalDataDir+'/'+'metabMap.csv', 'wb'))
    for key, value in namesDict.items():
        writer.writerow([key, value])
   
    return modelStatArray

################################################################################

# adjacencyListFromModel
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
    
################################################################################

# getModelStats
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
