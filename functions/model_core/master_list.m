
clear 

addpath(genpath('energyDeterminesMolecules-code'))

cd 'energyDeterminesMolecules-code'

% author: CORNELIUS BERGMANN
%
% modified 10.01.2024
%
% This script calculations mRNA and protein distributions with related
% costs for a list of parameter combinations and saves parameters and 
% energies in one large table.
% 
% This script uses the Parallel Computing Toolbox. To run the code without
% it, change the "parfor" in "run_list.m" to "for". Note that this will
% sharply increase the run time of the script.

%% set string to store data under  -  CHOOSE ONE

% no protein transport (trProb_p = 0)
% -----------------------------------------
% dataStr  = '2024_02_03-dendriteLength250'; denLength = 250;                
% dataStr  = '2024_02_03-dendriteLength500'; denLength = 500;    
% dataStr  = '2024_02_03-dendriteLength750'; denLength = 750;    
% dataStr  = '2024_02_03-dendriteLength1000'; denLength = 1000;    

% with protein transport (trProb_p = 0.1)
% -----------------------------------------
dataStr  = '2024_02_04-dendriteLength250'; denLength = 250;                
% dataStr  = '2024_02_04-dendriteLength500'; denLength = 500;    
% dataStr  = '2024_02_04-dendriteLength750'; denLength = 750;    
% dataStr  = '2024_02_04-dendriteLength1000'; denLength = 1000;    

%% parameter overview

% dendrite
%
% variable  |      values      | comments
% -------------------------------------------------------------------
% denLength |   250:250:1000   | 
% rhoSpines |         1        | 

denLength = denLength;                
rhoSpines = 1;

% transport model 
%
% variable |      values      | comments
% -------------------------------------------------------------------
% transp_m |     [0, 0.1]     | 
% transp_p |        0         | 
% v        |        1         | 
% beta     |        1         | 

v        = 1;
beta     = 1;
trProb_m = [0, 0.1]; 
trProb_p = 0.1; 

% mRNA and protein distribution
%
% variable |      values      | comments
% -------------------------------------------------------------------
% life_m   | [2, 8, 20] hours | 
% life_p   | [2, 8, 20] days  | 
% D_m      |  10.^(-4  :-2  ) | 
% D_p      |  10.^(-2.5:-0.5) | 
% eta_p    |  [10, 200, 5000] | 
% perm     |       1000       | 
% phi      |       0.95       | 
% tRate    |       0.01       | 
% somaRet  |      [0, 1]      | 

life_m   = [2, 8, 20].*3600; 
life_p   = [2, 8, 20].*3600.*24; 
D_m      = 10.^[-4, -3, -2];
D_p      = 10.^[-2.5, -1.5, -0.5];
eta_p    = [10, 200, 5000];      
permeab  = 1000;
phi      = 0.95;
tRate    = 0.01;
somaRet  = [0, 1];

% metabolic cost
%
% variable |      values       | comments
% -------------------------------------------------------------------
% n_aa     | [100, 500, 2000]  | 
% n_nt     | [3, 6, 15] x n_aa | 

n_aa     = [100,  500,  2000];   
n_nt     = [3, 6, 15];           

% bvp5c solver
%
% variable |      values       | comments
% -------------------------------------------------------------------
% maxRuns  |        10         | 

maxRuns  = 10;   

%% get and save ensemble mRNA diffusion constants 

% parameters used to fit mRNA diffusion constant (dendrite length and
% somatic concentration). Highly recommended to keep "length" at around 100
% and 100 is also a decent value for "atSoma"
length   = 100;
atSoma   = 100;
% Create matrix with one parameter combination per row. Relevant parameters
% are the instantaneous transport velocity, the switching rate out of the
% transported states, the mRNA half-life, the passive mRNA diffusion 
% constant and the fraction of transported mRNAs. For the latter, use only 
% positive values (zero equals no transport)
list     = {v, beta, trProb_m(trProb_m > 0), life_m, D_m};
tmp      = list;
[tmp{:}] = ndgrid(list{:});
list     = cell2mat(cellfun(@(m)m(:), tmp, 'uni', 0)); 
clear tmp 
% table for parameters and fitted diffusion constants
results  = table('Size'         , [size(list, 1), 6], ...
                 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'});
results.Properties.VariableNames = {'v [\mum/s]'         , ...
                                    'beta [1/s]'         , ...
                                    'transportedFraction', ...
                                    'half-life [s]'      , ...
                                    'D_0 [\mum^2/s]'     , ...
                                    'fittedD [\mum^2/s]' };

% iterate over mRNA diffusion constants and half-lives
for i = 1:size(list, 1)
    params            = [];
    % MUST BE: vPlus = vMinus
    params.vPlus      = list(i, 1);
    params.vMinus     = list(i, 1);
    % MUST BE: bPlus = bMinus
    params.bPlus      = list(i, 2);
    params.bMinus     = list(i, 2);
    params.transpProb = list(i, 3);
    params.alpha      = params.transpProb/(1 - params.transpProb)*(1/(1/params.bPlus + 1/params.bMinus));
    params.halflife   = list(i, 4);
    params.lambda     = log(2)/params.halflife;
    params.D          = list(i, 5);
    params.length     = length;
    params.atSoma     = atSoma;
    params.maxRuns    = maxRuns;
    % solve 3-state model and fit 1-state model
    out               = get_ensembleDiffConstFit(params);
    % store results in table
    results(i, :)     = {params.vPlus, ...
                         params.bPlus, ...
                         params.transpProb, ...
                         params.halflife, ...
                         params.D, ...
                         out.D_1State};
end
clear atSoma i length list out params 
% save ensemble mRNA diffusion constants
save(['files\', dataStr, '_fittedMRNADiffCoeffsTable.mat'], 'results')
clear results

%% get and save ensemble protein diffusion constants 

% parameters used to fit protein diffusion constant (dendrite length and
% somatic concentration). Highly recommended to keep "length" at around 100
% and 100 is also a decent value for "atSoma"
length   = 100;
atSoma   = 100;
% Create matrix with one parameter combination per row. Relevant parameters
% are the instantaneous transport velocity, the switching rate out of the
% transported states, the protein half-life, the passive protein diffusion 
% constant and the fraction of transported proteins. For the latter, use  
% only positive values (zero equals no transport)
list     = {v, beta, trProb_p(trProb_p > 0), life_p, D_p};
tmp      = list;
[tmp{:}] = ndgrid(list{:});
list     = cell2mat(cellfun(@(m)m(:), tmp, 'uni', 0)); 
clear tmp 
% table for parameters and fitted diffusion constants
results  = table('Size'         , [size(list, 1), 6], ...
                 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'});
results.Properties.VariableNames = {'v [\mum/s]'         , ...
                                    'beta [1/s]'         , ...
                                    'transportedFraction', ...
                                    'half-life [s]'      , ...
                                    'D_0 [\mum^2/s]'     , ...
                                    'fittedD [\mum^2/s]' };

% iterate over protein diffusion constants and half-lives
for i = 1:size(list, 1)
    params            = [];
    % MUST BE: vPlus = vMinus
    params.vPlus      = list(i, 1);
    params.vMinus     = list(i, 1);
    % MUST BE: bPlus = bMinus
    params.bPlus      = list(i, 2);
    params.bMinus     = list(i, 2);
    params.transpProb = list(i, 3);
    params.alpha      = params.transpProb/(1 - params.transpProb)*(1/(1/params.bPlus + 1/params.bMinus));
    params.halflife   = list(i, 4);
    params.lambda     = log(2)/params.halflife;
    params.D          = list(i, 5);
    params.length     = length;
    params.atSoma     = atSoma;
    params.maxRuns    = maxRuns;
    % solve 3-state model and fit 1-state model
    out               = get_ensembleDiffConstFit(params);
    % store results in table
    results(i, :)     = {params.vPlus, ...
                         params.bPlus, ...
                         params.transpProb, ...
                         params.halflife, ...
                         params.D, ...
                         out.D_1State};
end
clear atSoma  i length list out params 
% save ensemble protein diffusion constants
save(['files\', dataStr, '_fittedProteinDiffCoeffsTable.mat'], 'results')
clear results

%% create parameter list

% Create matrix with one parameter combination per row
list     = {trProb_m, trProb_p, somaRet, n_aa, n_nt, life_m, life_p, D_m, D_p, eta_p, denLength, rhoSpines, v, beta, permeab, phi, tRate, maxRuns};
tmp      = list;
[tmp{:}] = ndgrid(list{:});
list     = cell2mat(cellfun(@(m)m(:), tmp, 'uni', 0)); 
clear tmp trProb_m trProb_p somaRet n_aa n_nt life_m life_p D_m D_p 
clear eta_p denLength rhoSpines v beta permeab phi tRate maxRuns

%% remove undesired entries and save the list

% remove entries featuring active mRNA in combination with 100% somatic
% translation (it makes absolutely no sense to transport mRNAs if all of 
% them are retained in the soma)
removeInds          = list(:, 3) == 1 & list(:, 1) > 0;
list(removeInds, :) = [];
clear removeInds

% remove entries featuring passive mRNA in combination with > 0% dendritic
% translation (assuming that mRNAs in the dendrite effectively do not 
% leave the soma in the absence of active transport)
removeInds          = list(:, 3) < 1 & list(:, 1) == 0;
list(removeInds, :) = [];
clear removeInds

save(['files\', dataStr, '_parameterList.mat'], 'list')

%% define header for results table

header = {'transp_m'         ; ...
          'transp_p'         ; ...
          'somaRet'          ; ...
          'v'                ; ...
          'beta'             ; ...
          'D_m'              ; ...
          'D_p'              ; ...
          'ensembleD_m'      ; ...
          'ensembleD_p'      ; ...
          'tRate'            ; ...
          'halflife_m'       ; ...
          'halflife_p'       ; ...
          'permeability'     ; ...
          'length'           ; ...
          'rhoSpines'        ; ...
          'eta_p'            ; ...
          'phi'              ; ...
          'aminoAcids'       ; ...
          'nucleotides'      ; ...
          'maxRuns'          ; ...
          'mRNADendrTotal'   ; ...
          'mRNASomaTotal'    ; ...
          'mRNAInflux'       ; ...
          'proteinDendrTotal'; ...
          'proteinSpineTotal'; ...
          'proteinInflux'    ; ...
          'transportCost'    ; ...
          'translationCost'  ; ...
          'transcriptionCost'; ...
          'totalCost'        ; ...
          'totalCostInSoma'  ; ...
          'nMeshPts'         ; ...
          'nIterations'      ; ...
          'maxErr'           ; ...
          'totalCostPer10MicronSegment'};

%% get and save results 

[results, costPerSeg]               = run_list(size(list, 1), dataStr);
% transform to tables and add header
results                             = array2table(results);
results.Properties.VariableNames    = header(1:(end - 1));
costPerSeg                          = cell2table(costPerSeg);
costPerSeg.Properties.VariableNames = header(end);
% concatenate tables
results                             = [results, costPerSeg];
clear costPerSeg
save(['files\', dataStr, '_results.mat'], 'results')
clear header list
