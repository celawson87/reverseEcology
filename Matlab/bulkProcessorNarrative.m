%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bulkProcessor
% Copyright (c) 2014, Joshua J Hamilton and Katherine D McMahon
% Affiliation: Department of Bacteriology
%              University of Wisconsin-Madison, Madison, Wisconsin, USA
% URL: http://http://mcmahonlab.wisc.edu/
% All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each model should be in its own directory. Both xml and txt versions of
% the model are required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

% Retrieve the list of directories within the current directory
dataDir = 'ExternalData';
summaryDir = 'DataSummaries';
rawModelDir = 'RawModelFiles';
processedModelDir = 'ProcessedModelFiles';
if ~exist(strcat('../',processedModelDir), 'dir')
  mkdir(strcat('../',processedModelDir));
end

d = dir(strcat('../',rawModelDir));
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
% Perform post-processing
    modelStats = modelProcessorFunctionNarrative(strcat(nameFolds{i,1},'/',nameFolds{i,1}), dataDir, rawModelDir, processedModelDir);
% Store the model properties in the array
    results(i,:) = modelStats;
end

% Write results to file
fid = fopen(strcat('../',summaryDir,'/KBaseProcessing.txt'), 'w');
    fprintf(fid, 'Model\tGenes\tMetabs\tReactions\tBalanced\n');
for i = 1:dirSize(1)
    fprintf(fid, '%s\t%d\t%d\t%d\t%d\n', nameFolds{i,1}, results(i,1), results(i,2), results(i,3), results(i,4));
end
fclose(fid);