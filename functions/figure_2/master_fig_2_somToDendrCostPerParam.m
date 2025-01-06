
%% list of varied paramters for Figure 2 

varList      = cell2table({'D_m'                   , 'log_{10} D_m '                       , '[\mum^2/s]';   
                           'D_p'                   , 'log_{10} D_p '                       , '[\mum^2/s]';   
                           'halflife_m'            , 'mRNA half-life '                     , '[h]'       ;  
                           'halflife_p'            , 'Protein half-life '                  , '[d]'       ;  
                           'eta_p'                 , 'Protein copies/spine'                , ''          ; 
                           'aminoAcids'            , 'Amino acids/protein'                 , ''          ; 
                           'nucleotides:aminoAcids', {'# Non-coding/coding', 'nucleotides'}, ''          ; 
                           'length'                , 'Dendrite length'                     , '[\mum]'    });
% column heads
varList.Properties.VariableNames = {'Variable', 'Name for plots', 'Units for plots'};

%% load and prepare data

% date strings for the files with simulation results
dataStrList  = {'files\2024_02_03-dendriteLength250', ...
                'files\2024_02_03-dendriteLength500', ...
                'files\2024_02_03-dendriteLength750', ...
                'files\2024_02_03-dendriteLength1000'};
% consolidate simulation results
for indDate = 1:numel(dataStrList)
    % set date string
    dataStr  = dataStrList{indDate};
    % load data
    load([dataStr, '_results.mat'], 'results');
    % Consolidate and transform table to have one parameter combination per row
    datTmp   = get_consolidateParameterLists(results);
    clear results
    % remove cost per 50 microns, table concatenation with multiple
    % dendrite lengths does not allow it
    datTmp   = removevars(datTmp,{'totalCostPer10MicronSegment - som', 'totalCostPer10MicronSegment - den'});
    % Append data
    if indDate == 1
        data = datTmp;
    else
        data = [data; datTmp];
    end
end
clear dataStr datTmp indDate

%% plot cost ratio of somatic vs. dendritic mRNA

% take only one length
lengthToShow = 1000;
% choose which parameters' panels should be shown and in which order. 
varsToShow   = [3, 7, 5, 4, 6];

plot_somToDendrCostPerParam(data(data.('length') == lengthToShow, :), varList, varsToShow)
clear lengthToShow varsToShow
