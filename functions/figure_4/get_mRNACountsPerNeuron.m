function [sources, sourceDat] = get_mRNACountsPerNeuron()

% author: CORNELIUS BERGMANN
% 
% last modified 27.04.2023

%% Load data from Zeisel et al., 2015

% (https://doi.org/10.1126/science.aaa1934)
% load mRNA counts per gene and cell from the text file provided via 
% http://linnarssonlab.org/cortex
zeiselAll = readtable('Zeisel_2015.txt');
% keep only values from CA1 pyramidal cells. These cells are columns 683 to
% 1625, which we verified by the abundance of two marker genes (Spink8 and
% Fibcd1).
zeiselPyr = zeiselAll(:, 683:1625);
% compute average over cells
zeiselPyr = mean(zeiselPyr{:, 3:end}, 2, 'omitnan');
% remove NaNs and zeros
zeiselPyr = zeiselPyr(~isnan(zeiselPyr));
zeiselPyr = zeiselPyr(zeiselPyr > 0);
% correct for 20% sensitivity (Fig. S2)
zeiselPyr = 5.*zeiselPyr;

%% Load data from Perez et al., 2021

% (https://doi.org/10.7554/elife.63092)
% read sheets 'Dataset S3 - Dendrites' and 'Dataset S3 - Somata' 
[~, ~, perezDendr]                   = xlsread('Perez_2021.xlsx','Dataset S3 - Dendrites');
[~, ~, perezSoma]                    = xlsread('Perez_2021.xlsx','Dataset S3 - Somata');
% Load the gene names                                       in column "A"
%      the average mRNA counts in glutamatergic cells       in column "E"
perezDendr                           = cell2table(perezDendr(2:14844, [1, 5]));
perezSoma                            = cell2table(perezSoma( 2:14249, [1, 5]));
perezDendr.Properties.VariableNames  = {'gene name - Perez et al.', ...
                                        'mRNA count in dendrites - Perez et al.'};
perezSoma.Properties.VariableNames   = {'gene name - Perez et al.', ...
                                        'mRNA count in somata - Perez et al.'};
% array for matches between genes names of "Somata" and "Dendrites" tables
matches                              = nan(size(perezDendr, 1), 1);
% match gene names from somata and dendrites tables
for i = 1:size(perezDendr, 1)
    nameMatch          = strcmp(perezSoma.('gene name - Perez et al.'), perezDendr{i, 'gene name - Perez et al.'});
    matches(nameMatch) = perezSoma{nameMatch, 'mRNA count in somata - Perez et al.'};
end
% append matched soma counts to dendrite table
matches  = table(matches);
matches.Properties.VariableNames = {'mRNA count in somata - Perez et al.'};
perezGlu = [perezDendr, matches];
% add somatic and dendritic counts
perezGlu = perezGlu.('mRNA count in somata - Perez et al.') ...
         + perezGlu.('mRNA count in dendrites - Perez et al.');
% remove NaNs and zeros
perezGlu = perezGlu(~isnan(perezGlu));
perezGlu = perezGlu(perezGlu > 0);
% correct for 25% sensitivity (suppl. 1 to Fig. 1)
perezGlu = 4.*perezGlu;
clear i matches nameMatch perezDendr perezSoma

%% Set data sources and data

sources    = {'Zeisel et al., 2015', 'Perez et al., 2021'};
sourceDat  = {zeiselPyr, perezGlu};
