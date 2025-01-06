function [sources, sourceDat] = get_proteinCountsPerNeuron()

% author: CORNELIUS BERGMANN
% 
% last modified 27.04.2023

%% Load data from Helm et al., 2021

% (https://doi.org/10.1038/s41593-021-00874-w)
% read sheet 'final copy numbers per cell'
[~,~,helmProt]  = xlsread('Helm_2021.xlsx', 'Sheet1');
% Load the mean protein count per neuron in column C 
helmSpine       = helmProt(2:end, 3);
% discard 'non-detected' in row 48 and 'NaN' in the last two rows
helmSpine       = cellfun(@(x) strsplit(x, '±'), helmSpine([1:47, 49:(end-2)]), 'UniformOutput', false);
% get count per spine
helmMean        = zeros(size(helmSpine));
for i = 1:size(helmSpine, 1)
    tmp         = helmSpine{i, :};
    helmMean(i) = str2double(tmp{1});
end

%% Set data sources and data

sources    = {'Helm et al., 2021'};
sourceDat  = {helmMean};
