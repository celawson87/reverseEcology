%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finalModelCheck
% Copyright (c) 2014, Joshua J Hamilton and Katherine D McMahon
% Affiliation: Department of Bacteriology
%              University of Wisconsin-Madison, Madison, Wisconsin, USA
% URL: http://http://mcmahonlab.wisc.edu/
% All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each model should be in its own directory. Both xml and txt versions of
% the model are required.
% Code also requires:
%  table.tsv, a tab-delimited version of the SEED DB w/ subsystem info
%  ModelSEED-reactions-db.csv, a comma-separated version of the SEED DB w/
%  free energy infor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

% Retrieve the list of directories within the current directory
d = dir();
isub = [d(:).isdir];
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
dirSize = size(nameFolds);

% Create array to store results of PostProcessing
% Columns: genes, metabs, rxns, balanced (binary)
results = zeros(dirSize(1), 4);

% For each subdirectory ...
for i = 1:dirSize(1)
% Enter the subdirectory. Print the subdirectory name, # out of ##
    fprintf('Processing %s: %d of %d... \n', nameFolds{i,1}, i, dirSize(1));
% Perform a check of mass- and charge- balancing
    pathStr = (strcat(nameFolds{i,1},'/',nameFolds{i,1}, 'Balanced'));
    model = readCbModel(pathStr);
    
    fprintf('Checking mass- and charge-balancing ... \n');
    [massImbalance,imBalancedMass,imBalancedCharge,imBalancedBool] = checkMassChargeBalance(model);
    if any(imBalancedBool)
       fprintf('\nReactions remain unbalanced. \n');
        results(i,4) = 0;
    else
        fprintf('\nAll reactions are mass- and charge balanced\n');
        results(i,4) = 1;
    end
    
% Declare model parameters
    null = size(model.genes); results(i,1) = null(1);
    null = size(model.mets); results(i,2) = null(1);
    null = size(model.rxns); results(i,3) = null(1);
        
% Store the model properties in the array
end

% Write results to file
fid = fopen('Results.txt', 'w');
    fprintf(fid, 'Model\tGenes\tMetabs\tReactions\tBalanced\n');
for i = 1:dirSize(1)
    fprintf(fid, '%s\t%d\t%d\t%d\t%d\n', nameFolds{i,1}, results(i,1), results(i,2), results(i,3), results(i,4));
end
fclose(fid);