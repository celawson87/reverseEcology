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

%%% Import the subSystem data from SEED
fileID = fopen(strcat('../',dataDir,'/','table.tsv'));
mySeed = textscan(fileID, repmat('%q', 1, 7) , 'delimiter', '\t', 'CollectOutput', true);
fclose(fileID);
mySeed = mySeed{1,1};

%%% Import free energy data
fileID = fopen(strcat('../',dataDir,'/','ModelSEED-reactions-db.csv'));
dGData = textscan(fileID, repmat('%q', 1, 9) , 'delimiter', ',', 'CollectOutput', true);
fclose(fileID);
dGData = dGData{1,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the metabolite formulas from the text file
% We don't know how long a priori the metabolite list will be, so import
% the entire column and truncate.
fprintf('Importing charges ... \n');
model.metFormulas = importCharges(strcat('../',rawModelDir,'/',pathStr,'.txt'), 16, length(model.mets)+15);

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

% Import subsystem information
fprintf('Importing subsystem information ... \n');
    for j = 1:size(model.rxns,1)
        % Find the index in the Seed data matching the reaction
        myRxn = model.rxns{j};
        myRxn = regexprep(myRxn, '_[a-z]\d', '');
        myIndex = find(strcmp(mySeed(:,1), myRxn));
% Add subsystems
% If the entry exists, add information, otherwise write "None"        
        if ~isempty(myIndex)
% Check subsystem info
            mySub = mySeed(myIndex, 5);
            mySub = strrep(mySub,'"','');
            if regexp(mySub{1}, '[a-z]');
                model.subSystems(j) = mySub;
            else
                model.subSystems(j) = cellstr(['None']);
            end
% Check EC info
            myEC = mySeed(myIndex, 7);
            myEC = strrep(myEC,'"','');
            if regexp(myEC{1}, '[1-9]');
                model.rxnECNumbers(j) = myEC;
            else
                model.rxnECNumbers(j) = cellstr(['Unknown']);
            end
% Check if biomass        
        elseif regexp(myRxn, 'biomass')
            model.subSystems(j) = cellstr(['Biomass']);
            model.rxnECNumbers(j) = cellstr(['None']);
% Check if exchange
        elseif regexp(myRxn, 'EX_')
            model.subSystems(j) = cellstr(['Exchange']);
            model.rxnECNumbers(j) = cellstr(['None']);
        else
            model.subSystems(j) = cellstr(['Unknown']);
            model.rxnECNumbers(j) = cellstr(['Unknown']);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each reaction, look up the reaction in the 'THERMODYNAMIC FEASIBILTY'
% column. For each of the three cases, perform the appropriate action.

for j = 1:size(model.rxns,1)
    % Find the index in the Seed data matching the reaction
    myRxn = model.rxns{j};
    myRxn = regexprep(myRxn, '_[a-z]\d', '');
    myIndex = find(strcmp(dGData(:,1), myRxn));
    if myIndex
% If bidirectional, set model.rev = 1, set rxn bounds     
        if strcmp(dGData{myIndex, 9}, '<=>')
            model.rev(j) = 1;
            model.LB(j) = -1000;
            model.UB(j) = 1000;
% If forward-only, set model.rev = 0;
        elseif strcmp(dGData{myIndex, 9}, '=>')
            model.rev(j) = 0;
            model.LB(j) = 0;
            model.UB(j) = 1000;
% If reverse-only, switch reaction stoichiometry and set model.rev = 0;
        elseif strcmp(dGData{myIndex, 9}, '<=')
            model.S(:,j) = -model.S(:,j);
            model.rev(j) = 0;
            model.LB(j) = 0;
            model.UB(j) = 1000;
% If no entry, assume bidirectional, set model.rev = 1;        
        else 
            fprintf('No entry. Bidirectional by default.\n');
            model.rev(j) = 1;
            model.LB(j) = -1000;
            model.UB(j) = 1000;
        end
% If no entry, assume bidirectional, set model.rev = 1;     
    else
        model.rev(j) = 1;
        model.LB(j) = -1000;
        model.UB(j) = 1000;
    end
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