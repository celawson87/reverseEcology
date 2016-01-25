# From Genomes to Traits: Reverse Ecology of Uncultivated Freshwater Actinobacteria
### Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
### Department of Bacteriology
### University of Wisconsin-Madison, Madison, Wisconsin, USA
### [http://mcmahonlab.wisc.edu/](http://mcmahonlab.wisc.edu/)
### [https://github.com/joshamilton](https://github.com/joshamilton)
### All rights reserved.

***
## Data Folders

### refGenomes
Reference genomes used in this study. One folder per genome format.

### RawModelFiles
Contains raw output from KBase. One folder per genome. To construct a complete and valid SBML file from models built using the  Narrative interface, three files are required:

* AAA027E14.xml - SBML reconstruction from KBase  
* AAA027E14Compounds.tsv - Text format reconstruction from KBase  
* AAA027E14Reactions.tsv  - Text format reconstruction from KBase  

### excludedSamples
Raw and processed model files for genomes which have been excluded from this analysis.

### ExternalData
Additional files which are required to make various steps of the pipeline run.

* For reverseEcology.ipynb:
    * reactionsToRemove.txt - optional file. If you want to remove reactions from the models, list the reactions here.
    * taxonomy.csv - comma-separated file containing the following info for each sample: (ID, Lineage, Clade, Tribe, Color). Should be sorted alphabetically.
    * taxonomyTribes.csv - same as taxonomy.csv, but for all samples which are part of a tribe

* For pairwiseANI.ipynb:  
    * ANI_out - table of pairwise ANI calculations
    * COV_out - table of pairwise coverage calculations
    * taxonomySAG.csv - same as taxonomy.csv, for SAGs only

* For mergingGenomes.ipynb:  
    * taxonomySAG.csv - same as taxonomy.csv, for SAGs only
    * tribalColors.csv - sams as taxonomy.csv, but for tribes

* Also contains information about the genomes:
  * cogAbundance.csv - number of genes in each genome in each COG. Exported from IMG.
  * checkMMarkerGenes.csv -  presence/absence of marker genes used by CheckM
  * phylophlanMarkerGenes.csv - presence/absence of marker genes used by Phylophlan
  * phylosiftMarkerGenes.csv - presence/absence of marker genes used by Phylosift
  * masterColors.xlsx - spreadsheet showing coloration used in figures

### metadata
* 2015-08-21-FinalTree - Phylogenetic tree of genomes with final classifications
2015-08-18-actinoAlignment.phylip - Phylosift alignment used to generate the above tree
2015-08-18-actinoTreeCleanNames.nwk - RAxML-build tree file which I made pretty
* 2015-08-21-actinosMetadata(formatted).xlsx - metadata about the genomes, including estimated completeness, genome size, etc


***
## Code Folders

### Matlab
Codes for processing initial reconstructions from KBase. Deprecated in favor of Python codes.

### Perl
Codes for pushing genomes to KBase for model reconstruction. Deprecated since API access to KBase is broken (2016-01-25).

### Python
Codes for performing reverse ecology analysis.

* Function files called by workflow notebooks
  * metadataFunctions.py - functions for importing models and taxonomies
  * markerGeneFunctions.py - functions for computing completeness estimates from presence/absence of marker genes
  * sbmlFunctions.py - functions for working with SBML reconstructions from KBase  
  * graphFunctions.py - functions for working with graph representations of reconstructions  
  * seedFunctions.py - functions for computing reverse ecology metrics from seed sets
  * pairwiseANIFunctions.py - functions for performing ANI and coverage calculations
  * extraFunctions.py - extra functions which are currently unused in the final analysis

* Other Files
  * reverseEcology.py - master workflow for performing reverse ecology analysis  
  * mergingGraphs.py - master workflow for merging genome network graphs and performing reverse ecology analysis.

***
## Results Folders

### DataSummaries
Results of reverse ecology computations. Also contains a folder `MergedData` which contains results for analyses performed on tribe-level aggregates of genomes.

* ModelStatistics.txt - Numbers of genes, reactions, and metabolites in the SBML models.
* seedMatrixWeighted.csv - For each metabolite in the set of models, whether or not the compound is a seed (1 = yes). A metabolite may have value <1 if it can substitute for another compound.
* SeedCounts.txt - For ALL genomes, frequency with which a metabolite is a seed compound (1 = metabolite is a seed in 100% of genomes)
* MetabolicCompetitionScores.csv - metabolic competition between pairs of genomes
* MetabolicComplementarityScores.csv - metabolic complementarity between pairs of genomes
* metabolicCompetition.png - heatmap of competition scores
* metabolicComplementarity.png - heatmap of complementarity scores

### iPython Notebook
iPython notebook summaries of the project

* masterNotebook - master workflow for computing reverse ecology metrics.  
* ReverseEcologyOLD - old version of the above.

* genomesFromKBase.ipynb - workflow for converting SBML models from KBase into a form suitable for work with reverse ecology
* pairwiseANI - ANI- and coverage-based criteria for identifying tribes. Analysis of which partially-classified GFMs can be merged into tribes.  
* mergingGenomes - workflow for for merging network graphs of individual genomes into tribes.  
* samplingCompleteness.ipynb - accumulation curves of SAGS belonging to different tribes. Used to gauge how completely (or not) we have sampled the pan-genome of different tribes.
* toDo.ipynb - Wow, I have a lot of stuff to do!

### ProcessedModelFiles
Processed versions of KBase reconstructions.
* AAA023D18.xml - Reconstruction in SBML format
* AAA023D18RedAdjList.txt - Reconstruction in adjacency list format, used for RE calculations


Also includes results of seed set calculations:
* AAA023D18SeedSets.txt - list of seed compounds for the genome, one per line. If multiple seed compounds can trade-off, all are listed on a single line.
* AAA023D18SeedSetsWNames.txt - Same as above, but using metabolite names instead of IDs.
* AAA023D18SeedWeights.txt - list of all seed compounds with their weights
* AAA023D18SeedWeightsWNames.txt - same as above, but using metabolite names instead of IDs.

### MergedData
Same data as in ProcessedModelFiles folder, but for aggregates of genomes belonging to a single tribe.
