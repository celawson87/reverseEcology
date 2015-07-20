# From Genomes to Traits: Reverse Ecology of Uncultivated Freshwater Actinobacteria
##### Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
##### Department of Bacteriology
##### University of Wisconsin-Madison, Madison, Wisconsin, USA
##### http://mcmahonlab.wisc.edu/
##### All rights reserved.
***

#### Data Folders

##### ExternalData
Additional files which are required to make various steps of the pipeline run.

* For model pre-processing:
    * metaboliteTable.tsv - manual export of model SEED metabolite table
    * ModelSEED-reactions-db.csv - manual export of Model SEED database
    * reactionTable.tsv - manual export of model SEED reaction table

* For reverseEcology.ipynb:
    * reactionsToRemove.txt - optional file. If you want to remove reactions from the models, list the reactions here.
    * taxonomy.csv - comma-separated file containing the following info for each sample: (ID, Lineage, Clade, Tribe, Color). Should be sorted alphabetically.

* For pairwiseANI.ipynb:
    * ANI_out - table of pairwise ANI calculations
    * taxonomySAG.csv - same as taxonomy.csv, for SAGs only

* For mergingGenomes.ipynb:
    * taxonomySAG.csv - same as taxonomy.csv, for SAGs only
    * tribalColors.csv - sams as taxonomy.csv, but for tribes

* Auto-generated files
    * metabMap.csv - auto-generated file which maps SEED compound IDs to common names
    * reducedANI_OUT - ANI_out with information from SAGs only
    * betweenTribeANI - summary of max and min ANI of samples from different tribes
    * withinTribeANI - summary of max and min ANI of samples from same tribe

##### ProcessedModelFiles
Processed versions of reconstructions, including mass- and charge-balanced SBML files and graph representationsl.

##### RawModelFiles
Contains raw output from KBase. Each folder contains the reconstruction for a single genome. To construct a complete and valid SBML file from models built using the deprecated IRIS interface, three files are required:

* AAA023D18.xml - SBML reconstruction from KBase  
* AAA023D18.txt - Text format reconstruction from KBase  
* AAA023D18Annotation.txt - KBase annotation of the contigs.  

To construct a complete and valid SBML file from models built using the deprecated Narrative interface, three files are required:

* AAA027E14.xml - SBML reconstruction from KBase  
* AAA027E14Compounds.tsv - Text format reconstruction from KBase  
* AAA027E14Reactions.tsv  - Text format reconstruction from KBase  

##### excludedSamples
Raw and processed model files for genomes which have been excluded from this analysis.

##### DataSummaries
Summary statistics about genome-scale network reconstructions, metabolic network (di)graphs
and (di)graph representations. Results of reverse ecology computations, such as seed sets and competition scores. Contains a 'MergedData' subfolder containing the same information for merged reconstructions (e.g., of tribes). These data are not tracked, as they are re-generated each time the pipeline is run.

##### MergedData
Summary statistics about genome-scale network reconstructions, metabolic network (di)graphs
and (di)graph representations. These data are not tracked, as they are re-generated each time the pipeline is run.

#### Code Folders


##### iPython Notebook
iPython notebook summaries of the project

* ReverseEcology - master workflow for computing reverse ecology metrics.  
* pairwiseANI - ANI- and coverage-based criteria for identifying tribes. Analysis of which partially-classified GFMs can be merged into tribes.  
* mergingGenomes - workflow for for merging network graphs of individual genomes into tribes.  


##### Matlab
Codes for processing initial reconstructions from KBase.

* bulkProcessor.m - master function for processing models. Performs model-processing and simple mass-balancing. Creates an output file indicating which models require mass-balancing.  
* modelProcessorFunction.m - script which performs the actual model processing.  
* finalModelCheck.m - after models have been manually balanced, this script performs a final check of mass- and charge-balancing.  
* importCharges.m - COBRA Toolbox doesn't properly import charges from the SBML. This function imports them from the text version of the reconstruction.  
* reactionBalancer.m - simple script for checking mass- and charge-balancing on a single model.

##### Python
Codes for performing reverse ecology analysis.

* reverseEcology.py - master workflow for performing reverse ecology analysis  
* mergingGraphs.py - master workflow for merging genome network graphs and performing reverse ecology analysis.

* metadataFunctions.py - functions for importing models and taxonomies
* sbmlFunctions.py - functions for working with SBML reconstructions from KBase  
* graphFunctions.py - functions for working with graph representations of reconstructions  
* seedFunctions.py - functions for computing reverse ecology metrics from seed sets
* pairwiseANIFunctions.py - functions for performing ANI and coverage calculations
* extraFunctions.py - extra functions which are currently unused in the final analysis
