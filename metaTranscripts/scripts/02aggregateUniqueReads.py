###############################################################################
# countUniqueReads
# Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Process read count data from mapping of OMD-TOIL MT reads to our reference
# genomes.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

import HTSeq
import os

#%%#############################################################################
### Static folder structure
################################################################################
# Define fixed input and output files
genomeFolder = '../../refGenomes/concat'
sampleFolder = '../rawData'
mapFolder = '../bamfiles'
countFolder = '../htseq'
readsDir = countFolder+'/reads'

# Check that the new output directory exists and create if it doesn't
if not os.path.exists(countFolder):
        print "Creating output directory\n"
        os.makedirs(countFolder)

if not os.path.exists(readsDir):
        print "Creating output directory\n"
        os.makedirs(readsDir) 

#%%#############################################################################
### Read in MT and genome lists. Create DF to store read countDict.
################################################################################
# Read in list of MTs
mtList = []
for mt in os.listdir(sampleFolder):
    if mt.startswith('.'):
        next
    else:
       mtList.append(mt)

mtList = [mt.replace('.fastq', '') for mt in mtList]

# Read in list of genomes. Ignore internal standard genome.
genomeList = []
for genome in os.listdir(genomeFolder):
    if genome.endswith('.fna'):
       genomeList.append(genome)

genomeList = [genome.replace('.fna', '') for genome in genomeList]

#%%#############################################################################
### Run htseq-count
################################################################################

#for mt in mtList:
#    for genome in genomeList:
#        outFile = open(countFolder+'/'+mt+'-'+genome+'.CDS.out', "w")
#        subprocess.call(['/usr/local/bin/htseq-count',
#                         '-f', 'sam',
#                         '-r', 'pos',
#                         '-s', 'no',
#                         '-a', '0',
#                         '-t', 'CDS',
#                         '-i', 'locus_tag',
#                         '-m', 'intersection-strict',
#                         '-o', countFolder+'/'+mt+'-'+genome+'.CDS.sam',
#                         mapFolder+'/'+mt+'-'+genome+'.sam',
#                         genomeFolder+'/'+genome+'.gff'], stdout=outFile)
#        outFile.close()
#        
##%%#############################################################################
#### Obtain the same information for each (genome, COG) pairing
#################################################################################

# Portions of this code were borrowed verbatim from HTSeq-count
minQual = 0
featureType = 'CDS'
idAttr = 'locus_tag'
overlapMode = 'intersection-strict'

for mt in mtList:
    for genome in genomeList:
        samFile = mapFolder+'/'+mt+'-'+genome+'.sam'
        gffFile = genomeFolder+'/'+genome+'.gff'    

        # CReate feature array to store GFF data        
        features = HTSeq.GenomicArrayOfSets("auto", stranded = False)     
        
        # Create dictionary to store count data
        countDict = {}

        # Read in the GFF file
        gffObj = HTSeq.GFF_Reader(gffFile)   
        i = 0
        for gffLine in gffObj:
            if gffLine.type == featureType:
                featureID = gffLine.attr[idAttr]
                features[gffLine.iv] += featureID
            # Create dictionary entry
                countDict[gffLine.attr[idAttr]] = 0
            i = i + 1
        print(str(i)+' GFF lines processed')
            
        # Read in the SAM file
        readSeqFile = HTSeq.SAM_Reader(samFile)
        
        # Variables to store counts of read types
        alignOutCDS = 0
        nonUnique = 0
        i = 0   

        # Count reads
        for read in readSeqFile:
            i += 1
            # Skip ambiguous reads
            if read.optional_field("NH") > 1:
                nonUnique += 1
                continue
            iv_seq = (co.ref_iv for co in read.cigar if co.type == "M" and co.size > 0)
            fs = None
            for iv in iv_seq:
                for iv2, fs2 in features[iv].steps():
                    if len(fs2) > 0 or overlapMode == "intersection-strict":
                        if fs is None:
                            fs = fs2.copy()
                        else:
                            fs = fs.intersection(fs2)
            if fs is None or len(fs) == 0:
                alignOutCDS += 1
            else:
                countDict[list(fs)[0]] += 1
                # Write the read to file
                with open(readsDir+'/'+list(fs)[0]+".reads", "a") as myfile:
                   myfile.write(read.original_sam_line+"\n")
            
        print(str(i)+' SAM alignments processed')      
            
        with open(countFolder+'/'+mt+'-'+genome+'.CDS.out', "w") as outFile:
            for fn in sorted(countDict.keys()):   
                outFile.write("%s\t%d\n" % (fn, countDict[fn]))
            outFile.write("__align_outside_CDS\t%d\n" % alignOutCDS)
            outFile.write("__alignment_not_unique\t%d" % nonUnique)