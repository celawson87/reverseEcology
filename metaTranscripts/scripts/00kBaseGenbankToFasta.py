###############################################################################
# CodeTitle.py
# Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Code description.
################################################################################

#%%#############################################################################
### Import packages
################################################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os

#%%#############################################################################
### Define folder structure
################################################################################
gbkDir = '../../refGenomes/gbk'
fnaDir = '../../refGenomes/fna_kbase'

if not os.path.exists(fnaDir):
        print "Creating output directory\n"
        os.makedirs(fnaDir)
        
#%%#############################################################################
### Create list of genomes to process
################################################################################

genomeList = []
for genome in os.listdir(gbkDir):
    if genome.endswith('.gbk'):
        genomeList.append(genome)

genomeList = [genome.replace('.gbk', '') for genome in genomeList]

#%%#############################################################################
### Convert Genbank files to FASTA nucleotide files
### These are necessary to run clustering algorithms which require input gene
### sequences.
################################################################################

for genome in genomeList:
    inFile = open(gbkDir+'/'+genome+'.gbk')
    outFile = open(fnaDir+'/'+genome+'.fna', "w")
    
    for gbkRecord in SeqIO.parse(inFile, "genbank"):
        print >>outFile, ">"+gbkRecord.id
        print >>outFile, gbkRecord.seq
 
    inFile.close()
    outFile.close()
