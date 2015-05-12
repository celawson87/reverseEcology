###############################################################################
# extraFunctions.py
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Extra code snippets used for analyses which didn't make it to the final
# analysis.
################################################################################


################################################################################

#%%
# Identify the most highly-connected metabolites
myPct = 0.01
removeMetabs = gf.findTopMetab(myPct, graphList)

#%%
# Remove such metabolites from the graphs of each network, and rerun graphStatistics
reducedGraphStatArray = np.empty([numSubDir, 4], dtype = int)

count = 0
reducedGraphFile = open('ReducedGraphStatistics.txt', 'w')
reducedGraphFile.write('Model,Nodes,Edges,Total Components,Size of Largest\n')

reducedGraphList = []
reducedDiGraphList = []
for curDir in dirList:
    print 'Processing graph', count+1, 'of', numSubDir, ':', curDir
# Read in adjacency list and convert to graph object
    myGraph = graphList[count]
    myGraph.remove_nodes_from(removeMetabs)
    myDiGraph = diGraphList[count]
    myDiGraph.remove_nodes_from(removeMetabs)
    
    reducedGraphList.append(myGraph)
    reducedDiGraphList.append(myDiGraph)
# Read graph statistics                       
    reducedGraphStatArray[count:] = gf.getGraphStats(myGraph)
    reducedGraphFile.write('%s,%i,%i,%i,%i\n' % (curDir, 
                                                 reducedGraphStatArray[count,0], 
                                       reducedGraphStatArray[count,1], 
                                       reducedGraphStatArray[count, 2],
                                       reducedGraphStatArray[count, 3] ) )
    count = count + 1
reducedGraphFile.close()

################################################################################
