###############################################################################
# mtFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Set of functions for working with metatranscriptomes.
################################################################################ 

# Import Python packages.
import cobra
import networkx as nx
import pandas as pd
import re

# These custom-written modules should have been included with the package
# distribution. 
import metadataFunctions as mf

################################################################################ 

# seedsToCOGs
# A function which identifies COGs associated with metabolism of seed
# compounds. Inputs are the following directories:

# dirList: A list of models to perform this analysis on
# genomeModelDir: Directory containing network and sbml models for individual
# genomes
# mergedModelDir: Directory contaiing netowrk models and seed lists for 
# consensus genome models
# taxonFile: File giving taxonomic information for individual genomes

# Output is a file genomeSeedsToCogs.txt which lists all COGs associated with
# metabolism of each seed compound. One entry per line.

def compoundsToCOGs(modelList, genomeModelDir, mergedModelDir, taxonFile, metabType):

    # Create a list of genomes belonging to each linage/clade/tribe
    lineageDict =  mf.importTaxonomy(taxonFile, 'Lineage')
    cladeDict =  mf.importTaxonomy(taxonFile, 'Clade')
    tribeDict =  mf.importTaxonomy(taxonFile, 'Tribe')
    sampleDict = dict(lineageDict.items() + cladeDict.items() + tribeDict.items())
    
    # Loop over each merged model...
    
    for curDir in modelList:
    # Read in its seed compounds from the weighted list. Drop the weight and 
    # add a column for the COGs.
    
        seedDF = pd.read_csv('../'+mergedModelDir+'/'+curDir+'/'+curDir+metabType+'.txt', header=None, names=['Metab'], index_col=0)
    #    seedDF = seedDF.drop('1.000000', 1)
        seedDF = pd.concat([seedDF,pd.DataFrame(columns=['COGs'])])
    
        seedToCogDict = dict.fromkeys(seedDF.index, [])
    
    # Read its adjacency list representation
        mergedDiGraph = nx.read_adjlist('../'+mergedModelDir+'/'+curDir+'/'+curDir+'AdjList.txt',
                                create_using=nx.DiGraph())    
    
    # Retrieve the list of genomes associated with the merged model
        genomeList = sampleDict[curDir]
    
    # For each genome ... 
        for genome in genomeList:
    
        # Read in its edge mappings
            genomeEdgeMapping = {}
            with open('../'+genomeModelDir+'/'+genome+'/'+genome+'RxnEdges.txt', mode ='r') as inFile:
                for line in inFile:
                    edgeMap = line.strip().split('\t')
                    genomeEdgeMapping[edgeMap[0]+','+edgeMap[1]] = edgeMap[2]
    
        # Read in the SBML model containing GPRs
            genomeModel = cobra.io.read_sbml_model('../'+genomeModelDir+'/'+genome+'/'+genome+'.xml')
        
        # Read in the COG dictionary
            cogDict = {}
            inFile = open('../'+genomeModelDir+'/'+genome+'/'+genome+'COGs.txt', 'r')
            for line in inFile.readlines():
                key = line.strip().split(',')[0]
            # Reformat the key for compatability with COBRA
            # Replace .genome. with .
                key = re.sub('\.genome\.', '.', key)
            # Replace . with underscore
                key = re.sub('\.','_', key)
                value = line.strip().split(',')[1]
                cogDict[key] = value
            inFile.close()
        
        # For each seed compound in the consensus genome, find all sets of outward-pointing arcs
        # For each such, look up all associated reactions in the genome's adjacency list
        # For all reactions, look up all the associated genes in the model
        # For all genes, look up the proper COG
            for seed in seedToCogDict.keys():
            # Look up all sink nodes
                sinkNodeList = mergedDiGraph[seed].keys()
    
            # Generate the list of associated reactions
                reactionList = []
                for node in sinkNodeList:
                    if seed+','+node in genomeEdgeMapping.keys():
                        reactionList.append(genomeEdgeMapping[seed+','+node])
            # Reduce to unique entries
                reactionList = list(set(reactionList))
            
            # Generate the list of associated genes
                geneList = []
                for reaction in reactionList:
                    for gene in genomeModel.reactions.get_by_id(reaction).genes:
                        if gene.id != 'Unknown':                    
                            geneList.append(gene.id)
            # Reduce to unique entries
                geneList = list(set(geneList))
            
                cogList = []
                for gene in geneList:
                    cogList.append(cogDict[gene])
            # Reduce to unique entries
                cogList = list(set(cogList))
            
            # Update seedToCogDict with new COG elements
                seedToCogDict[seed] = seedToCogDict[seed] + cogList
    
    # Now that the dict has been fully populated, reduce keys to lists of unique values
        for seed in seedToCogDict.keys():
            seedToCogDict[seed] = list(set(seedToCogDict[seed]))
        
    # Write the results to file
        outFile = open('../'+mergedModelDir+'/'+curDir+'/'+curDir+metabType+'ToCogs.txt', mode ='w')
        for key, value in seedToCogDict.items():
            outFile.write(key+':'+','.join(value)+'\n')
        outFile.close()

    return