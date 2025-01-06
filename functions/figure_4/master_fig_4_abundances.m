
%% list of varied paramters for Figures 4A,B,C

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

%% load and prepare data for Figures 4A,B,C
     
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

%% Load experimental data for Figure 4B

% derive mRNA counts per neuron from published datasets
[sources_m, sourceDat_m] = get_mRNACountsPerNeuron();
% multiply data from Zeisel et al., 2015 with 20, leading to a total
% scaling (compared to raw data) of 100 (multiplication with 5 was
% performed to account for 20% sensitivity)
isZeisel              = strcmp(sources_m, 'Zeisel et al., 2015');
zeiselDat             = sourceDat_m{isZeisel};
zeiselDat             = 20.*zeiselDat;
sourceDat_m{isZeisel} = zeiselDat;
clear isZeisel zeiselDat

% derive protein counts per neuron from published datasets
[sources_p, sourceDat_p] = get_proteinCountsPerNeuron();

%% Figures 4A & 4B

% limits of plots for panels A, B and boxplot insets
yLim.A = [-1, 19];
yLim.B = [-1, 13];
boxLim = [-1, 14];
plot_abundancesCheapestStrategy_withExpData(data, yLim, boxLim, sources_m, sourceDat_m, sources_p, sourceDat_p)
clear boxLim yLim

%% Figure 4C

% take only one length
lengthToShow = 250;
% plot mRNA and protein abundance, clustered along the cheapest strategy
plot_abundancesDendriticVsSomatic(data, lengthToShow)
clear lengthToShow

%% Load experimental data for Figure 4D

% load mRNA abundance for somata-enriched and neurite-enriched genes 
[sources_m, sourceDatSom_m, sourceDatNeur_m, signif_m] = get_mRNAAbundanceVsNeuriteEnrichment();
% derive data from Zappulo et al., 2017
isZappulo       = strcmp(sources_m, 'Zappulo et al., 2017');
zappuloSom_m    = sourceDatSom_m{isZappulo};
zappuloNeur_m   = sourceDatNeur_m{isZappulo};
zappuloSignif_m = signif_m{isZappulo};
clear isZappulo sources_m sourceDatNeur_m sourceDatSom_m signif_m

% load protein abundance for somata-enriched and neurite-enriched genes 
[sources_p, sourceDatSom_p, sourceDatNeur_p, signif_p] = get_proteinAbundanceVsNeuriteEnrichment();
% derive data from Zappulo et al., 2017
isZappulo       = strcmp(sources_p, 'Zappulo et al., 2017');
zappuloSom_p    = sourceDatSom_p{isZappulo};
zappuloNeur_p   = sourceDatNeur_p{isZappulo};
zappuloSignif_p = signif_p{isZappulo};
clear isZappulo sources_p sourceDatNeur_p sourceDatSom_p signif_p

%% Figure 4D

yLim_m = [-0.75, 5.25]; 
yLim_p = [6.25, 11.75]; 
plot_experimentalAbundancesVsNeuriteEnrichment(zappuloSom_m, zappuloNeur_m, zappuloSignif_m, yLim_m, ...
                                               zappuloSom_p, zappuloNeur_p, zappuloSignif_p, yLim_p)
clear yLim_m yLim_p
