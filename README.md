# Auxotrophies and Interactomes in Mixed Culture of Freshwater Bacteria
### Copyright (c) 2016, Joshua J Hamilton and Katherine D McMahon
### Department of Bacteriology
### University of Wisconsin-Madison, Madison, Wisconsin, USA
### [http://mcmahonlab.wisc.edu/](http://mcmahonlab.wisc.edu/)
### [https://github.com/joshamilton](https://github.com/joshamilton)
### All rights reserved.

***
## Data Folders

### refGenomes
Reference genomes used in this study. One folder per genome format.

### sbml

#### sbml/raw

Contains raw output from KBase. One folder per genome. To construct a complete and valid SBML file from models built using the  Narrative interface, three files are required:

* acIB2-FNE-F8.xml - SBML reconstruction from KBase  
* acIB2-FNE-F8Compounds.tsv - Text format reconstruction from KBase  
* acIB2-FNE-F8Reactions.tsv  - Text format reconstruction from KBase  

#### sbml/processed

Processed versions of KBase reconstructions.
* acIB2-FNE-F8.xml - Reconstruction in SBML format

***
## Code Folders

### python
Codes for performing reverse ecology analysis.

* metadataFunctions.py - functions for importing models and taxonomies
* sbmlFunctions.py - functions for working with SBML reconstructions from KBase  
* graphFunctions.py - functions for working with graph representations of reconstructions  
* seedFunctions.py - functions for computing reverse ecology metrics from seed sets

### scripts
Scripts for working with KBase.
* loadGenomes - script for pushing your refGenomes to KBase
* concatGbk.sh - concatenate gbk files for each contig into a single file for each genome
* kBaseGenbankToFasta.py - convert gbk to ffn and ffa format
* kBaseGenbankToGff.py - convert gbk to gff format
* cleanUpGFF.pl - remove extra comments from gff files

***
## Results Folders

### data
Results of reverse ecology computations.

* modelStats.tsv - Numbers of genes, reactions, and metabolites in the SBML models.
* metabMap.csv - A mapping between SEED compound IDs and metabolite names
* MetabolicCompetitionScores.csv - metabolic competition between pairs of genomes
* MetabolicComplementarityScores.csv - metabolic complementarity between pairs of genomes
* metabolicCompetition.png - heatmap of competition scores
* metabolicComplementarity.png - heatmap of complementarity scores

Also contains folders with seed set calculations for  individual genomes. Each folder includes:
  * acIB2-FNE-F8AdjList.txt - Reconstruction in adjacency list format
  * acIB2-FNE-F8RedAdjList.txt - Reduced reconstruction in adjacency list format, used for RE calculations
  * acIB2-FNE-F8SeedSets.txt - list of seed compounds for the genome, one per line. If multiple seed compounds can trade-off, all are listed on a single line.
  * acIB2-FNE-F8SeedSetsWNames.txt - Same as above, but using metabolite names instead of IDs.
  * acIB2-FNE-F8SeedWeights.txt - list of all seed compounds with their weights
  * acIB2-FNE-F8SeedWeightsWNames.txt - same as above, but using metabolite names instead of IDs.

### jupyter
Jupyter notebook summaries of the project

* FNE-F8 - master workflow for analysis of the FNE-F8 mixed culture.
