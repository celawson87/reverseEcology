###############################################################################
# cobraModelProcessor.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
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
################################################################################

# Import Python modules 
import cobra
import fileinput
import numpy as np
import pandas as pd
import re
import os

# Import Python modules 
# These custom-written modules should have been included with the package
# distribution. 
import metadataFunctions as mf

# Define local folder structure for data input and processing.
processedDataDir = 'ProcessedModelFiles'
rawModelDir = 'RawModelFiles'
summaryStatsDir = 'DataSummaries'

# Check that folders exist and create them if necessary
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

    ### Read in external data
    # Import metabolite charges
    cpdData = pd.read_csv('../'+rawModelDir+'/'+curDir+'/'+curDir+'Compounds.tsv', delimiter='\t', index_col=0)

################################################################################                   

    # Before reading in the SBML file, update gene loci so they read in properly
    # KBase generates gene loci of the form kb|g.######.CDS.###
    # Transform them to be of the form curDir_CDS_###
    for myLine in fileinput.FileInput('../'+rawModelDir+'/'+curDir+'/'+curDir+'.xml', inplace = True):
        print re.sub('kb\|g\.\d+\.CDS\.(\d+)', curDir+'_CDS_\g<1>', myLine).strip()
    fileinput.close()

################################################################################                   
    
    # Read in model from SBML
    model = cobra.io.read_sbml_model('../'+rawModelDir+'/'+curDir+'/'+curDir+'.xml')

################################################################################                   

    # Remove undesired reactions, including:
    badRxnList = []
    print 'Removing bad reactions'
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
    model.remove_reactions(badRxnList, delete=True, remove_orphans=True)                        

################################################################################                   

    # Update the metabolite formulas
    print 'Updating metabolite formulas'
    for curMetab in model.metabolites:
        # Retrieve the metabolite name w/o compartment info
        # Look up the appropriate index in the cpdData DF
        curMetab.formula = cpdData.loc[re.sub('_[a-z]\d', '', curMetab.id)][1]

################################################################################                   

    # Check mass- and charge- balancing   
    imbalCounter = 0
    print 'Correcting mass- and charge-balancing'

    # Check for reactions known to be imbalanced and manually correct them
    # If metabolite cpd03422 exists, update its charge to +1
    for curMetab in model.metabolites:
        if curMetab.id == 'cpd03422_c0':
            curMetab.charge = 1

    for curRxn in model.reactions:
    # If reaction rxn07295 exists, update its stoichiometry
        if curRxn.id == 'rxn07295_c0':
            curRxn.reaction = 'cpd00007_c0 + cpd00033_c0 <=> cpd00025_c0 + 3.0 cpd00067_c0 + cpd14545_c0'
    # If reaction rxn08808 exists, update its stoichiometry        
        elif curRxn.id == 'rxn08808_c0':
            curRxn.reaction = 'cpd00001_c0 + cpd15341_c0 <=> cpd00067_c0 + cpd00908_c0 + cpd01080_c0'

    # Heuristic for proton balancing
    for curRxn in model.reactions:
        imbalDict = curRxn.check_mass_balance()
        if len(imbalDict) != 0:
            # If imbalancing due to protons alone, correct it
            if imbalDict['H'] == imbalDict['charge']:
                curRxn.add_metabolites({model.metabolites.get_by_id('cpd00067_c0'): -1*imbalDict['H']})
            # Otherwise, indicate unbalanced reactions remain
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
    
    print 'Writing to file'
    cobra.io.write_sbml_model(model, '../'+processedDataDir+'/'+curDir+'/'+curDir+'.xml')
    
# Write the results to file
modelSizeDF.to_csv('../'+summaryStatsDir+'/modelStats.tsv', sep='\t')