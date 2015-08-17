function [results] = modelProcessorFunction(pathStr, dataDir, rawModelDir, processedModelDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AddFormulas.m
% Copyright (c) 2014, Joshua J Hamilton and Katherine D McMahon
% Affiliation: Department of Bacteriology
%              University of Wisconsin-Madison, Madison, Wisconsin, USA
% URL: http://http://mcmahonlab.wisc.edu/
% All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draft reconstructions from Kbase require some post-processing. This
% script does several important things:
% 1. Import metabolite charges
% 2. Check mass- and charge-balancing of reactions in the reconstruction
% 3. Import subsystem information
% 4. Import free energy information. Calculate free energy range and
% reverse reaction direction if appropriate (to comply w/ COBRA convention
% of reaction being forward-only or reversible).
% 5. Remove biomass, exchange, spontaneous, DNA/RNA biosynthesis reactions
% 6. Remove the corresponding genes
% 7. Modify metabolite names to comply with COBRA naming conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function returns an array of four integers:
% Columns: genes, metabs, rxns, balanced (binary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Establish vector for results
results = zeros(4,1);
% Read in model from SBML
fprintf('Reading in the model ... \n');
model = readCbModel(strcat('../',rawModelDir,'/',pathStr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the metabolite formulas from the text file
% Read in the TSV file w/ formula and charge info
cpdData = tdfread(strcat('../',rawModelDir,'/',pathStr,'Compounds.tsv'));

% Look up model metabs in the TSV file and update formula field:
model.metFormulas = cell(length(model.mets),1);
for i=1:length(model.mets)
% Regex to retrieve compound without compartment info
    myCmpd = model.mets{i};
    myCmpd = regexprep(myCmpd, '_[a-z]\d', '');
% Look up index in cpdData.name
    index = strmatch(myCmpd, cpdData.id);
% Update metFormulas
    myFormula = cpdData.formula(index,:);
    myFormula = strtrim(myFormula);
    model.metFormulas{i} = myFormula; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check mass- and charge-balancing
fprintf('Checking mass- and charge-balancing ... \n');
[massImbalance,imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model);
if any(imBalancedBool)
    fprintf('\nReactions are unbalanced. Attempting to balance on protons alone.\n');
    % Perform balancing
    H_index = find(ismember(model.mets,'cpd00067_c0'));
        for i = 1:length(model.rxns)
        model.S(H_index, i) = model.S(H_index, i) - full(massImbalance(i,1));
        end
        
    % Perform manual balancing of reactions known to be unbalanced in KBase
    % If metabolite cpd03422 exists, update its charge to +1
    cpdIndex = find(not(cellfun('isempty', strfind(model.mets, 'cpd03422'))));
    if not(isempty(cpdIndex))
        model.metCharge(cpdIndex) = 1;
    end
   
    % If reaction rxn07295 exists, update its stoichiometry
    rxnIndex = find(not(cellfun('isempty', strfind(model.rxns, 'rxn07295'))));
    if not(isempty(rxnIndex))
        model.S(findMetIDs(model,'cpd00007_c0'), rxnIndex) = -1;
        model.S(findMetIDs(model,'cpd00025_c0'), rxnIndex) = 1;
        model.S(findMetIDs(model,'cpd00067_c0'), rxnIndex) = 1;
    end

    % If reaction rxn08808 exists, update its stoichiometry
    rxnIndex = find(not(cellfun('isempty', strfind(model.rxns, 'rxn08808'))));
    if not(isempty(rxnIndex))
        model.S(findMetIDs(model,'cpd00067_c0'), rxnIndex) = 1;
        model.S(findMetIDs(model,'cpd12547_c0'), rxnIndex) = 0;
        model.S(findMetIDs(model,'cpd15341_c0'), rxnIndex) = -1;
    end
    
   % Check mass- and charge-balancing again
    [massImbalance,imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model);
    if any(imBalancedBool)
       fprintf('\nReactions remain unbalanced. \n');
       results(4) = 0;
    else
        fprintf('\nAll reactions are mass- and charge balanced\n');
        results(4) = 1;
    end 
else
    fprintf('\nAll reactions are mass- and charge balanced\n');
    results(4) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Eliminating reactions without genes ... \n');
% First identify reactions without GPRs (exchange reactions, BM reaction)
emptyCells = cellfun('isempty', model.grRules);
emptyRxns = model.rxns(emptyCells);
model = removeRxns(model, emptyRxns);

% Then identify reactions whose GPR is wholly "unknown"
emptyCells = find(strcmp('Unknown', model.grRules));
emptyRxns = model.rxns(emptyCells);
model = removeRxns(model, emptyRxns);

fprintf('Removing protein, DNA, and RNA biosynthesis ... \n');
removeList = cellstr(['rxn13782_c0'; 'rxn13783_c0'; 'rxn13784_c0']);
model = removeRxns(model, removeList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write to file
fprintf('Establishing COBRA compliance ... \n');
% Make Kbase model consistent with expecations of COBRA Toolbox
model.mets = regexprep(model.mets, '_([a-z])\d', '[$1]');
% Remove trailing 0s
model.rxns = regexprep(model.rxns, 'biomass0', 'biomass');
model.rxns = regexprep(model.rxns, '_[a-z]\d', '');
% rxnNames
model.rxnNames = regexprep(model.rxnNames, '\s[a-z]\d', '');
% metNames
model.metNames = regexprep(model.metNames, '_[a-z]\d', '');

% pathStr contains both the folder and file name. Split along '/' to
% identify the relevant dir.
splitStr = regexp(pathStr, '/', 'split');
dir = splitStr(1);
if ~exist(strcat('../',processedModelDir,'/',dir{1,1}), 'dir')
  mkdir(strcat('../',processedModelDir,'/',dir{1,1}));
end

writeCbModel(model, 'sbml', strcat('../',processedModelDir,'/',pathStr,'Balanced.xml'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in again, to ensure unused genes are removed
fprintf('Eliminating unused genes ... \n');
model = readCbModel(strcat('../',processedModelDir,'/',pathStr,'Balanced.xml'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare model parameters
null = size(model.genes); results(1,1) = null(1);
null = size(model.mets); results(2,1) = null(1);
null = size(model.rxns); results(3,1) = null(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform final write to file
fprintf('Performing final write ... \n');

writeCbModel(model, 'sbml', strcat('../',processedModelDir,'/',pathStr,'Balanced.xml'));