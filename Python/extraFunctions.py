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

def findTopMetab(myPct, graphList):

    aggNodeCount = Counter()
    for graph in graphList:
        aggNodeCount = aggNodeCount + Counter(nx.degree(graph))
    
    aggNodeList = aggNodeCount.most_common()
    # Indices for splitting metabolites into two sets.
    totalNodes = len(aggNodeCount)
    pctIndex = int(math.ceil(totalNodes*myPct))
    
    # Plot the number of edges associated with each metabolite
    x = np.linspace(1, totalNodes, totalNodes)
    y = zip(*aggNodeList)[1]
    x0 = x[0:pctIndex-1]
    x1 = x[pctIndex:totalNodes]
    y0 = y[0:pctIndex-1]
    y1 = y[pctIndex:totalNodes]
    
    pyplot.loglog(x0, y0, marker='.', color='red', linestyle='none')
    pyplot.loglog(x1, y1, marker='.', color='black', linestyle='none')
    pyplot.xlim(0, len(x))
    pyplot.ylim(0, max(y))
    pyplot.xlabel('Metabolite Rank')
    pyplot.ylabel('Number of Edges')
    
    # Return a list of the most-connected metabolites
    metabFile = open('connectedMetabs.txt', 'w')
    metabFile.write('Metabolite, Edges\n')
    
    removeMetabs = []
    for count in range(pctIndex):
        metabFile.write('%s,%i\n' % (aggNodeList[count][0], 
                                     aggNodeList[count][1] ) )
        removeMetabs.append(aggNodeList[count][0])
                                 
    metabFile.close()

    return removeMetabs
