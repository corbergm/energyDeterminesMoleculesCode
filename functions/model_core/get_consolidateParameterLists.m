function datOut = get_consolidateParameterLists(datIn)

% author: CORNELIUS BERGMANN
%
% modified 03.02.2024
%
% This function accepts raw simulation results and consolidates the output into a
% more useful table format.
% A priori, each row corresponds to one simulation with a given localization
% pathway. Hence, each parameter combination appears in multiple rows, once 
% for each localization pathway. This function transforms this into a table
% with one parameter combination and all corresponding results (with all
% localization pathways) in one row.

% Before:
% parameters 1 | somatic   translation             | results (cost etc.)
% parameters 1 | dendritic translation             | results (cost etc.)

% After:
% parameters 1 | results with dendr. transl. | results with som. transl. 

%% split data into one table per localization pathway

% column with nucleotides:amino acids
datTmp  = table(datIn.('nucleotides')./datIn.('aminoAcids'));
datTmp.Properties.VariableNames = {'nucleotides:aminoAcids'};
% append to data
datIn   = [datIn, datTmp];
clear datTmp
% split tables along localization pathways
datDen  = datIn(datIn.('somaRet') == 0 & datIn.('transp_m') == 0.1, :);
datSom  = datIn(datIn.('somaRet') == 1                            , :);
clear datIn
% sort tables to get the same order of rows within each of them
datDen  = sortrows(datDen ,{'D_m', 'D_p', 'halflife_m', 'halflife_p', 'eta_p', 'aminoAcids', 'nucleotides:aminoAcids', 'length'}, ...
                           {'descend', 'descend', 'descend', 'descend', 'descend', 'descend', 'descend', 'descend'});
datSom  = sortrows(datSom ,{'D_m', 'D_p', 'halflife_m', 'halflife_p', 'eta_p', 'aminoAcids', 'nucleotides:aminoAcids', 'length'}, ...
                           {'descend', 'descend', 'descend', 'descend', 'descend', 'descend', 'descend', 'descend'});

%% merge tables 

% derive columns with parameters which are similar for all localization
% pathways
paramTab = datDen(:, ["v", "beta", "tRate", "permeability", "rhoSpines", "phi", "maxRuns", ...
                      "D_m", "D_p", "halflife_m", "halflife_p", "eta_p", "aminoAcids", "nucleotides", "nucleotides:aminoAcids", "length"]);
% derive columns with parameters which are different among localization
% pathways
datDen   = datDen(:, ["transp_m", "transp_p", "somaRet", "ensembleD_m", "ensembleD_p", ...
                      "mRNADendrTotal", "mRNASomaTotal", "mRNAInflux", "proteinDendrTotal", "proteinSpineTotal", "proteinInflux", ...
                      "transportCost", "translationCost", "transcriptionCost", "totalCost", "totalCostInSoma", "totalCostPer10MicronSegment", ...
                      "nMeshPts", "nIterations", "maxErr"]);
datSom   = datSom(:, ["transp_m", "transp_p", "somaRet", "ensembleD_m", "ensembleD_p", ...
                      "mRNADendrTotal", "mRNASomaTotal", "mRNAInflux", "proteinDendrTotal", "proteinSpineTotal", "proteinInflux", ...
                      "transportCost", "translationCost", "transcriptionCost", "totalCost", "totalCostInSoma", "totalCostPer10MicronSegment", ...
                      "nMeshPts", "nIterations", "maxErr"]);
% add localization pathway names to pathway-specific variable names
for varInd = 1:size(datDen, 2)
    datDen.Properties.VariableNames{varInd} = [datDen.Properties.VariableNames{varInd}, ' - den' ];
    datSom.Properties.VariableNames{varInd} = [datSom.Properties.VariableNames{varInd}, ' - som'    ];
end
clear varInd
% concatenate parameters similar and different among localization pathways
datOut     = [paramTab, datSom, datDen];
clear paramTab datSom datDen
