
% This script fits ensemble mRNA and protein diffusion constants using
% a range of transport model, mRNA/protein parameters, and dendrite
% lengths to validate if using 100 microns for the fitting procedure is a
% valid choice. To this end, several checks of the diffusion constant fits
% are shown, including the relative change of ensemble diffusion constant
% with dendrite length, the solver error, the number of mesh points used by
% the solver, the ratio of anterograde to retrograde transport, the actual
% vs required ratio of transported molecules, and, finally, the error of
% the least-squares fitting of ensemble diffusion constants.

% set date string  
dataStr = '2024_02_01_transportModel'; 

%% generic parameter overview

% transport model parameters
%
% variable |      values                 | comments
% -------------------------------------------------------------------
% v        |      [0.2, 0.5, 2, 5]       | reference value: 1
% beta     |      [0.2, 0.5, 2, 5]       | reference value: 1
% transp_m | [0.02, 0.05, 0.2, 0.3, 0.5] | reference value: 0.1
% length   |            100              | 

v           = 1;                           % [micron/s]; 
vList       = [0.2; 0.5; 2; 5];            % [micron/s]; 
beta        = 1;                           % [1/s]; 
betaList    = [0.2; 0.5; 2; 5];            % [1/s]; 
trProb      = 0.1;                         % [%]; 
trProbList  = [0.02; 0.05; 0.2; 0.3; 0.5]; % [%]; 
length      = 100;                         % [micron]; 

% mRNA parameters 
%
% variable |      values      | comments
% -------------------------------------------------------------------
% life_m   | [2, 8, 20] hours | 
% D_m      |    10.^(-4:-1)   | in this setting we also consider 0.1

life        = [2, 8, 20].*3600;           % [s]; 
D           = [0.0001, 0.001, 0.01, 0.1]; % [micron^2/s]; 

%% table with the reference values vor v, beta, transport probability

% "reference" transport model parameter values
refVals = [v, beta, trProb];
% header for "generic" values
header  = {'v'                    ; ...
           'beta'                 ; ...
           'transport probability'};
refVals = array2table(refVals);
refVals.Properties.VariableNames = header;
clear header

%% create parameter list

% create a matrix with all combinations of the transport model parameters.
% Each parameter is varied independently
transpList = [v                        , beta                        , trProb                      ;
              vList                    , beta.*ones(size(vList))     , trProb.*ones(size(vList))   ;
              v.*ones(size(betaList))  , betaList                    , trProb.*ones(size(betaList));
              v.*ones(size(trProbList)), beta.*ones(size(trProbList)), trProbList                  ];
clear beta betaList v vList trProb trProbList

% Create matrix with all combinations of diffusion constant, half-life, and
% dendrite length
mRNAList           = {life, D, length};
tmp               = mRNAList;
[tmp{:}]          = ndgrid(mRNAList{:});
mRNAList           = cell2mat(cellfun(@(m)m(:), tmp, 'uni', 0)); 
clear tmp D length life

% create a list with all combinations of transport model parameters and
% mRNA parameters
list = [];
for ind = 1:size(mRNAList, 1)
    list = [list; [transpList, repmat(mRNAList(ind, :), size(transpList))]];
end
clear ind mRNAList transpList
% save list and the transport model reference values
save(['files\', dataStr, '_parameterList.mat'], 'list', 'refVals')

%% prepare for parfor usage

% slice list into single columns for improved "parfor" performance
vList      = list(:, 1);
betaList   = list(:, 2);
trProbList = list(:, 3);
lifeList   = list(:, 4);
DList      = list(:, 5);
lenList    = list(:, 6);
% number of parameter combinations
nCmbs      = size(list, 1);
% initialise matrix for results with 14 columns per entry
results    = zeros(nCmbs, 14);
clear list

%% define header for results table

header = {'v'                    , ...
          'beta'                 , ...
          'transport probability', ...
          'halflife'             , ...
          'D'                    , ...
          'length'               , ...
          'atSoma'               , ...
          'maxRuns'              , ...
          'D_1State'             , ...
          'err_3State'           , ...
          'nMeshPts_3State'      , ...
          'anteroToRetro'        , ...
          'anteroRetroToTotal'   , ...
          'err_fitting'          };

%% fit diffusion constants for parameter list

parfor i = 1:nCmbs
    % suppress solver detail output
    warning('off', 'MATLAB:bvp5c:RelTolNotMet');
    % show progress
    disp(i)
    % set parameters
    params            = [];
    params.length     = lenList(i);
    % MUST BE: vPlus = vMinus
    params.vPlus      = vList(i);
    params.vMinus     = vList(i);
    % MUST BE: bPlus = bMinus
    params.bPlus      = betaList(i);
    params.bMinus     = betaList(i);
    params.halflife   = lifeList(i);
    params.lambda     = log(2)/params.halflife;
    params.D          = DList(i);
    params.transpProb = trProbList(i);
    params.alpha      = params.transpProb/(1 - params.transpProb)*(1/(1/params.bPlus + 1/params.bMinus));
    params.atSoma     = 100;
    params.maxRuns    = 15;
    out               = get_ensembleDiffConstFit(params);
    % store results
    results(i, :)     = [params.vPlus                , ...
                         params.bPlus                , ...
                         params.transpProb           , ...
                         params.halflife             , ...
                         params.D                    , ...
                         params.length               , ...
                         params.atSoma               , ...
                         params.maxRuns              , ...
                         out.D_1State                , ...
                         out.err_3State              , ...
                         out.nMeshPts_3State         , ...
                         mean(out.anteroToRetro)     , ...
                         mean(out.anteroRetroToTotal), ...
                         out.err_fitting             ]
end
clear betaList lenList vList bList trProbList lifeList DList nCmbs out

%% transform to table and add header

results                          = array2table(results);
results.Properties.VariableNames = header;
% save results
save(['files\', dataStr, '_results.mat'], 'results')
clear refVals header results
