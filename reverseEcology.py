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

# Import additional functions which are found in relevant scripts
import sbmlFunctions

# Retrieve listing of model subdirectories
dirList = filter(os.path.isdir, os.listdir(os.getcwd()))
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
    modelStatArray[count:] = sbmlFunctions.getModelStats(model)
    sbmlFunctions.adjacencyListFromModel(model)
    count = count + 1

# Write model statistics to file
colLabels = ['Genes', 'Metabolites', 'Reactions']
modelStatDF = pandas.DataFrame(modelStatArray, columns=colLabels, index=dirList) 
modelStatDF.to_csv('ModelStatistics.txt')
