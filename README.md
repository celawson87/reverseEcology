# From Genomes to Traits: Reverse Ecology of Uncultivated Freshwater Actinobacteria
##### Copyright (c) 2015, Joshua J Hamilton and Katherine D McMahon
##### Department of Bacteriology
##### University of Wisconsin-Madison, Madison, Wisconsin, USA
##### http://http://mcmahonlab.wisc.edu/
##### All rights reserved.
***

#### Data Folders

DataSummaries - Summary statistics about genome-scale models and (di)graph representations used in this project

ExternalData - SBML files from KBase are missing important information. Additional files are downloaded from the Model SEED to provide this information.

* ModelSEED-reactions-db.csv - Information about model reactions. "Thermodynamic Feasibility" is used to determine proper reaction direction based on dG values.  
* table.tsv - Information about model reactions. "Subsystem" column is used to add additional information to the SBML file.

ProcessedModelFiles - Processed versions of reconstrucitons, including mass- and charge-balanced SBML files and graph representations.

RawModelFiles - Contains raw output from KBase. Each folder contains the reconstruction for a single genome. To construct a complete and valid SBML file, three files are required:

* AAA023D18.xml - SBML reconstruction from KBase  
* AAA023D18.txt - Text format reconstruction from KBase  
* AAA023D18Annotation.txt - KBase annotation of the contigs.  


#### Code Folders


iPython Notebook - iPython notebook summary of the project

Matlab - Codes for processing initial reconstructions from KBase.

* bulkProcessor.m - master function for processing models. Performs model-processing and simple mass-balancing. Creates an output file indicating which models require mass-balancing.  
* modelProcessorFunction.m - script which performs the actual model processing.  
* finalModelCheck.m - after models have been manually balanced, this script performs a final check of mass- and charge-balancing.  
* importCharges.m - COBRA Toolbox doesn't properly import charges from the SBML. This function imports them from the text version of the reconstruction.  
* reactionBalancer.m - simple script for checking mass- and charge-balancing on a single model.

Python - Codes for performing reverse ecology analysis.

* reverseEcology.py - master function for performing reverse ecology analysis  
* sbmlFunctions.py - functions for working with SBML reconstructions from KBase  
* graphFunctions.py - functions for working with graph representations of reconstructions  
* extraFunctions.py - extra functions which are currently unused in the final analysis

