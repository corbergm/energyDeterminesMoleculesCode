function [results, costPerSeg] = run_list(nCmbs, dataStr)

% author: CORNELIUS BERGMANN
%
% modified 03.02.2024
%
% This function solves the model equation system (detailed in
% "get_distribution") for the list of parameters from the file
% "[dataStr, '_parameterList.mat']". 
% Ensemble mRNA and protein diffusion constants are loaded from
% "[dataStr, '_fittedMRNADiffCoeffsTable.mat']" and
% "[dataStr, '_fittedProteinDiffCoeffsTable.mat']"
% 
% Simulation results are stored in the "results" matrix 
% Input: 
% nCmbs      - number of parameter combinations to be solved
% dataStr    - data string to identify parameter list and obtained results
%
% Output: 
% results      - nCmbs x 34 matrix containing simulation results, for column
%                details see below
% costPerSeg   - nCmbs x 1 cell array containing the total cost per 10 micron
%                dendrite segement

% initialise matrix for results with 34 columns per entry
results = zeros(nCmbs, 34);
% initialise cell array for the total cost per 10 micron dendrite segement
costPerSeg = cell(nCmbs, 1);

parfor j = 1:nCmbs
    
    % suppress solver detail output
    warning('off', 'MATLAB:bvp5c:RelTolNotMet');

    % show progress
    disp(j)
    % load parameter list
    tmp                 = load(['files\', dataStr, '_parameterList.mat'], 'list');
    list                = tmp.list;
    % create parameter structure
    params              = struct();
    params.transp_m     = list(j,  1);
    params.transp_p     = list(j,  2);
    params.somaRet      = list(j,  3);
    params.aminoAcids   = list(j,  4);   
    params.nucleotides  = list(j,  4)*list(j, 5);
    params.halflife_m   = list(j,  6);    
    params.halflife_p   = list(j,  7);
    params.D_m          = list(j,  8);
    params.D_p          = list(j,  9);
    params.eta_p        = list(j, 10);
    params.length       = list(j, 11);    
    params.rhoSpines    = list(j, 12); 
    params.v            = list(j, 13); 
    params.beta         = list(j, 14); 
    params.permeability = list(j, 15); 
    params.phi          = list(j, 16); 
    params.tRate        = list(j, 17); 
    params.maxRuns      = list(j, 18); 

    %% get ensemble mRNA diffusion coefficient
    
    if params.transp_m > 0
        % load table with baseline ("passive") diffusion coefficients, mRNA 
        % halflife and fitted ensemble diffusion coefficients
        fitted_D_m    = load(['files\', dataStr, '_fittedMRNADiffCoeffsTable.mat'], 'results');
        fitted_D_m    = fitted_D_m.results;
        % get the correct ensemble diffusion coefficient from "fitted_D_m"
        matchV        = fitted_D_m.('v [\mum/s]')          == params.v;
        matchBeta     = fitted_D_m.('beta [1/s]')          == params.beta;
        matchTransp_m = fitted_D_m.('transportedFraction') == params.transp_m;
        matchD_m      = fitted_D_m.('D_0 [\mum^2/s]')      == params.D_m;
        matchLife_m   = fitted_D_m.('half-life [s]')       == params.halflife_m;
        fitted_D_m    = fitted_D_m.('fittedD [\mum^2/s]');
        params.DEns_m = fitted_D_m(matchV & matchBeta & matchTransp_m & matchD_m & matchLife_m);
    else
        params.DEns_m = params.D_m;
    end
    
    %% get ensemble protein diffusion coefficient
    
    if params.transp_p > 0
        % load table with baseline ("passive") diffusion coefficients, protein 
        % halflife and fitted ensemble diffusion coefficients
        fitted_D_p    = load(['files\', dataStr, '_fittedProteinDiffCoeffsTable.mat'], 'results');
        fitted_D_p    = fitted_D_p.results;
        % get the correct ensemble diffusion coefficient from "fitted_D_p"
        matchV        = fitted_D_p.('v [\mum/s]')           == params.v;
        matchBeta     = fitted_D_p.('beta [1/s]')           == params.beta;
        matchTransp_p = fitted_D_p.('transportedFraction') == params.transp_p;
        matchD_p      = fitted_D_p.('D_0 [\mum^2/s]')       == params.D_p;
        matchLife_p   = fitted_D_p.('half-life [s]')        == params.halflife_p;
        fitted_D_p    = fitted_D_p.('fittedD [\mum^2/s]');
        params.DEns_p = fitted_D_p(matchV & matchBeta & matchTransp_p & matchD_p & matchLife_p);
    else
        params.DEns_p = params.D_p;
    end
    
    %% get mRNA and protein distributions 
    
    solution = get_distribution(params);

    %% Derive mRNA and protein distribution features  
    
    % initialise structures for mRNA and protein
    mRNA                = [];
    protein             = [];
    % Set the points at which the distributions should be evaluated,
    % distance should be 0.5 microns.
    xval                = linspace(0, params.length, 2*params.length);
    % Derive the mRNA and protein distributions. Round to 10 digits to 
    % delete errors below mass conservation accuracy that disturb the plots.
    mRNA.dendr          = round(deval(solution, xval, 2), 10);
    protein.dendr       = round(deval(solution, xval, 4), 10);
    % derive the protein distribution in spines along the dendrite
    protein.spine       = (protein.dendr.*params.permeability)./((protein.dendr.*params.permeability)./(params.rhoSpines*params.eta_p) + 1);
    % integrate dendritic mRNA and dendritic protein
    mRNA.dendrTotal     = trapz(xval, mRNA.dendr);
    protein.dendrTotal  = trapz(xval, protein.dendr); 
    % integrate spine protein
    protein.spineTotal  = trapz(xval, protein.spine); 
    % Derive the mRNA and protein influx
    [~, mRNA.influx]    = deval(solution, 0, 2);
    [~, protein.influx] = deval(solution, 0, 4);
    mRNA.influx         = - params.DEns_m*mRNA.influx;
    protein.influx      = - params.DEns_p*protein.influx;
    % derive absolute somatic mRNAs
    mRNA.somaTotal      = protein.influx/params.tRate;

    %% Get metabolic cost

    cost = get_cost(params, mRNA, protein, xval);
    
    %% Save results
    
    costPerSeg{j} = cost.totalPerSeg;

    results(j, :) = [params.transp_m     , ...
                     params.transp_p     , ...
                     params.somaRet      , ...
                     params.v            , ...
                     params.beta         , ...
                     params.D_m          , ...
                     params.D_p          , ...
                     params.DEns_m       , ...
                     params.DEns_p       , ...
                     params.tRate        , ...
                     params.halflife_m   , ...
                     params.halflife_p   , ...
                     params.permeability , ...
                     params.length       , ...
                     params.rhoSpines    , ...
                     params.eta_p        , ...
                     params.phi          , ...
                     params.aminoAcids   , ...
                     params.nucleotides  , ...
                     params.maxRuns      , ...
                     mRNA.dendrTotal     , ...
                     mRNA.somaTotal      , ...
                     mRNA.influx         , ...
                     protein.dendrTotal  , ...
                     protein.spineTotal  , ...
                     protein.influx      , ...
                     cost.transp         , ...
                     cost.transl         , ...
                     cost.transcr        , ...
                     cost.total          , ...
                     cost.totalInSoma    , ...
                     solution.nMeshPts   , ...
                     solution.nIter      , ...
                     solution.maxErr     ];
end
