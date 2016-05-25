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
import cobra.core.Formula
import csv
import fileinput
import numpy as np
import os
import pandas as pd
import re

# Import custom Python modules 
import metadataFunctions as mf

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
        model = cobra.io.read_sbml_model('../'+processedDataDir+'/'+curDir+'/'+curDir+'.xml')

# Create dictionary of metabolite names
        for metab in model.metabolites:
            namesDict[metab.id] = metab.name

# Update description field
        model.id = curDir;

# Read model statistics by invoking sbmlFunctions.getModelStats
        modelStatArray[count:] = getModelStats(model)
        modelFile.write('%s,%i,%i,%i\n' % ('../'+processedDataDir+'/'+curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )

# Create adjacency list and write to file
        adjacencyListFromModel(model, processedDataDir)
        reactionEdgesFromModel(model, processedDataDir)
        count = count + 1

# Close files containing summary data
    modelFile.close()

# Write completed dictionary to file as a csv file ExternalData/metabMap.csv
    writer = csv.writer(open('../'+externalDataDir+'/'+'metabMap.csv', 'wb'))
    for key, value in namesDict.items():
        writer.writerow([key, value])
   
    return modelStatArray

################################################################################

# dirListToAdjacencyListWithRemoval
# This function iterates over a list of genome directories and converts each
# genome scale model from an SBML file to an adjacency list. Adjacency lists 
# for each genome-scale model are written as text files in each genome 
# directory. Summary statistics about each graph are written in the
# summaryStatsDir as well. Removes from each model the reaction specified in
# externalDataDir/reactionsToRemove.txt

def dirListToAdjacencyListWithRemoval(dirList, externalDataDir, processedDataDir, summaryStatsDir):

    numSubDir = len(dirList)

# Read in the list of reactions to remove
    with open('../'+externalDataDir+'/'+'reactionsToRemove.txt') as rxnFile:
        rxnList = [x.strip() for x in rxnFile.readlines()]
    
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
        model = cobra.io.read_sbml_model('../'+processedDataDir+'/'+curDir+'/'+curDir+'.xml')

# Remove reactions slated for removal
        for rxn in rxnList:
# Check if the reaction exists in the model. If so, delete it.
            if rxn in model.reactions:
                model.remove_reactions(rxn, delete=True, remove_orphans=True)
            
# Create dictionary of metabolite names
        for metab in model.metabolites:
            namesDict[metab.id] = metab.name

# Update description field
        model.id = curDir;

# Read model statistics by invoking sbmlFunctions.getModelStats
        modelStatArray[count:] = getModelStats(model)
        modelFile.write('%s,%i,%i,%i\n' % ('../'+processedDataDir+'/'+curDir, modelStatArray[count,0], 
                                    modelStatArray[count,1], 
                                    modelStatArray[count, 2] ) )

# Create adjacency list and write to file
        adjacencyListFromModel(model, processedDataDir)
        reactionEdgesFromModel(model, processedDataDir)
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
    myFile = open('../'+processedDataDir+'/'+model.id+'/'+model.id+'AdjList.txt', 'w')

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

# reactionEdgesFromModel
# Function to convert a cobrapy model object to a list of (source, sink) pairs
# for reaction in the model. For each reaciton, creates an edge between all 
# (reactant, product) pairs and indicates the reaction. If a reaction is 
# reversible, also creates edges between all (product, reactant) pairs.
# Input: cobrapy model object, model directory
# Output: None.

def reactionEdgesFromModel(model, processedDataDir):

# Establish a file for the adjacency list
    myFile = open('../'+processedDataDir+'/'+model.id+'/'+model.id+'RxnEdges.txt', 'w')

# For each reaction, loop over the reactants. For each reactant, loop over the 
# reaction products and create an edge between the reactant and products. If a 
# reaction is reversible, repeat the process in reverse, creating an edge
# between each product and reactant. Also record reaction associated with each
# edge.
    for myRxn in model.reactions:
        for myReactant in myRxn.reactants:
            for myProduct in myRxn.products:
                myFile.write(myReactant.id+'\t')
                myFile.write(myProduct.id+'\t')
                myFile.write(myRxn.id+'\n')
        if myRxn.reversibility == True:
            for myProduct in myRxn.products:
                for myReactant in myRxn.reactants:
                    myFile.write(myProduct.id+'\t')
                    myFile.write(myReactant.id+'\t')
                    myFile.write(myRxn.id+'\n')
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

################################################################################

# Draft reconstructions from Kbase require some post-processing. This script 
# does several important things:
# 1. Reformat gene locus tags
# 2. Remove biomass, exchange, spontaneous, DNA/RNA biosynthesis reactions and 
# their corresponding genes
# 3. Import metabolite formulas
# 4. Check mass- and charge-balancing of reactions in the reconstruction
# 5. Remove trailing 0s from reaction and metabolite names

# The post-processing has a major shortcoming. When KBase detects that one or 
# more subunits of a complex are present, it creates a "full" GPR by adding 
# 'Unknown' genes for the other subunits. CobraPy currently lacks functions to 
# remove the genes. As such, these model should not be used to perform any 
# simulations which rely on GPRs.

# Each model should be in its own directory in the 'RawModelFiles' folder. Both
# SBMl and TSV versions from KBase are required.

# As output, the code returns processed SBML files in the 'processedDataDir'
# folder. Also returns a summary of the model sizes, in the 'summaryStatsDir'
# folder.

def processSBMLforRE(rawModelDir, processedDataDir, summaryStatsDir):

# # Check that folders exist and create them if necessary
    if not os.path.exists('../'+processedDataDir):
        os.makedirs('../'+processedDataDir)
    if not os.path.exists('../'+summaryStatsDir):
        os.makedirs('../'+summaryStatsDir)
    
# Import the list of models
    dirList = mf.getDirList('../'+rawModelDir)
    numSubDir = len(dirList)

# Create an array to store results
# Columns: genes, metabs, rxns, balanced (binary)
    modelSizeDF = pd.DataFrame(index = dirList, columns=['Genes', 'Metabolites', 'Reactions', 'Balanced'])

# Intialize a counter
    count = 1
    for curDir in dirList:
    # Establish vector for results
        resultArray = np.zeros((4,1), dtype=int);

    # Print the subdirectory name
        print 'Processing model ' + curDir + ', ' + str(count) + ' of ' + str(numSubDir)

    # Import metabolite charges
        cpdData = pd.read_csv('../'+rawModelDir+'/'+curDir+'/'+curDir+'Compounds.tsv', delimiter='\t', index_col=0)

################################################################################                   

    # Before reading in the SBML file, update gene loci so they read in properly
    # In the old models ...
    # KBase generates gene loci of the form kb|g.######.CDS.###
    # Transform them to be of the form curDir_CDS_###

    # In the new models ...
    # KBase generates gene loci of the form curDir.genome.CDS.###
    # Transform them to be of the form curDir_CDS_###
    
        for myLine in fileinput.FileInput('../'+rawModelDir+'/'+curDir+'/'+curDir+'.xml', inplace = True):
#            print re.sub('kb\|g\.\d+\.CDS\.(\d+)', curDir+'_CDS_\g<1>', myLine).strip()
            print re.sub(curDir+'\.genome\.CDS\.(\d+)', curDir+'_CDS_\g<1>', myLine).strip()
        fileinput.close()

################################################################################                   
    
    # Read in model from SBML
        model = cobra.io.read_sbml_model('../'+rawModelDir+'/'+curDir+'/'+curDir+'.xml')

################################################################################                   

    # Remove undesired reactions, including:
        badRxnList = []
#        print 'Removing bad reactions'
        for curRxn in model.reactions:
        # Exchange reactions
            if re.match('EX_', curRxn.id):
                badRxnList.append(curRxn)
        # Reactions w/o GPRs
            elif len(curRxn.gene_reaction_rule) == 0:
                badRxnList.append(curRxn)
        # Protein, DNA, and RNA synthesis
            elif curRxn.id == 'rxn13782_c0' or curRxn.id == 'rxn13783_c0' or curRxn.id == 'rxn13784_c0':
                badRxnList.append(curRxn)
        # Spontaneous reactions, whose GPR is fully 'unknown'
            elif curRxn.gene_reaction_rule == 'Unknown':
                badRxnList.append(curRxn)     
        # Transport reactions, based on keywords
            elif re.search('transport', curRxn.name) or re.search('permease', curRxn.name) or re.search('symport', curRxn.name) or re.search('diffusion', curRxn.name) or re.search('excretion', curRxn.name) or re.search('export', curRxn.name) or re.search('secretion', curRxn.name) or re.search('uptake', curRxn.name) or re.search('antiport', curRxn.name):
                badRxnList.append(curRxn)
        # Transport reactions which don't get picked up based on keywords
            elif curRxn.id == 'rxn05226_c0' or curRxn.id == 'rxn05292_c0' or curRxn.id == 'rxn05305_c0' or curRxn.id == 'rxn05312_c0' or curRxn.id == 'rxn05315_c0' or curRxn.id == 'rxn10945_c0':
                badRxnList.append(curRxn)
        model.remove_reactions(badRxnList, delete=True, remove_orphans=True)                        

        print 'The remaining extracellular metabolites are:'
        for curMetab in model.metabolites:
            if re.search('_e0', curMetab.id):
                print curMetab.id

################################################################################                   

    # Update the metabolite formulas
#        print 'Updating metabolite formulas'
        for curMetab in model.metabolites:
        # Retrieve the metabolite name w/o compartment info
        # Look up the appropriate index in the cpdData DF
            curMetab.formula = cobra.core.Formula.Formula(cpdData.loc[re.sub('_[a-z]\d', '', curMetab.id)][1])
            curMetab.formula.id = cpdData.loc[re.sub('_[a-z]\d', '', curMetab.id)][1]

################################################################################                   

    # Check mass- and charge- balancing   
        imbalCounter = 0
#        print 'Correcting mass- and charge-balancing'

    # Check for reactions known to be imbalanced and manually correct them
    # If metabolite cpd03422 exists, update its charge to +1
        for curMetab in model.metabolites:
            if curMetab.id == 'cpd03422_c0':
                curMetab.charge = 1

        for curRxn in model.reactions:
    # If reaction rxn07295 exists, update its stoichiometry
            if curRxn.id == 'rxn07295_c0':
                curRxn.reaction = 'cpd00007_c0 + cpd00033_c0 <=> cpd00025_c0 + 3.0 cpd00067_c0 + cpd14545_c0'
                print 'Manually correcting an imbalance'
    # If reaction rxn08808 exists, update its stoichiometry        
            elif curRxn.id == 'rxn08808_c0':
                curRxn.reaction = 'cpd00001_c0 + cpd15341_c0 <=> cpd00067_c0 + cpd00908_c0 + cpd01080_c0'
                print 'Manually correcting an imbalance'
    # If reaction rxn12822 exists, update its stoichiometry        
            elif curRxn.id == 'rxn12822_c0':
                curRxn.reaction = '2.0 cpd00023_c0 + cpd11621_c0 <=> cpd00024_c0 + cpd00053_c0 + 2.0 cpd00067_c0 + cpd11620_c0'
                print 'Manually correcting an imbalance'

    # Heuristic for proton balancing
        for curRxn in model.reactions:
            imbalDict = curRxn.check_mass_balance()
            if len(imbalDict) != 0:
            # If imbalancing due to protons alone, correct it
                if imbalDict['H'] == imbalDict['charge']:
                    curRxn.add_metabolites({model.metabolites.get_by_id('cpd00067_c0'): -1*imbalDict['H']})
                    print 'Re-balancing on the basis of protons'
                else:
                    imbalCounter = imbalCounter + 1                
                    print 'Reaction ' + str(curRxn.id) + ' remains unbalanced'                
                
    # Inform of results
        if imbalCounter != 0:
            modelSizeDF.loc[curDir][3] = 0
        else:
            modelSizeDF.loc[curDir][3] = 1
            print 'All reactions are balanced'              
        
################################################################################            
    
    ### Update names to remove trailing zeros
        for curComp in model.compartments:
            model.compartments[curComp] = re.sub('_\d', '', model.compartments[curComp])
            model.compartments[re.sub('\d', '', curComp)] = model.compartments.pop(curComp)
        
        for curMetab in model.metabolites:
            curMetab.id = re.sub('\d$', '', curMetab.id)
            curMetab.name = re.sub('_[a-z]\d$', '', curMetab.name)
            curMetab.compartment = re.sub('\d$', '', curMetab.compartment)
        
        for curRxn in model.reactions:
            curRxn.id = re.sub('\d$', '', curRxn.id)
            curRxn.name = re.sub('_[a-z]\d$', '', curRxn.name)
                
################################################################################                   

    # Store the model properties in the array, write the model output, and increase the counter
        modelSizeDF.loc[curDir][0] = len(model.genes)
        modelSizeDF.loc[curDir][1] = len(model.metabolites)
        modelSizeDF.loc[curDir][2] = len(model.reactions)

    # Perform final write to file
    # Check that output dir exists and create if necessary    
        if not os.path.exists('../'+processedDataDir+'/'+curDir):
            os.makedirs('../'+processedDataDir+'/'+curDir)
    
#        print 'Writing to file'
        cobra.io.write_sbml_model(model, '../'+processedDataDir+'/'+curDir+'/'+curDir+'.xml')
        count = count + 1
    
# Write the results to file
    modelSizeDF.to_csv('../'+summaryStatsDir+'/modelStats.tsv', sep='\t')

    return
