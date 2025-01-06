

% author: CORNELIUS BERGMANN
% 
% last modified 12.12.2022

% This script calculates mRNA and protein distributions and metabolic cost
% for two localization strategies for a parameter set inspired by CaMKIIa 
% mRNA and protein and finally plots the metabolic costs.

% suppress solver detail output
warning('off', 'MATLAB:bvp5c:RelTolNotMet');

%% CaMKIIa-inspired parameters for mRNA and protein

% generic parameters
params              = [];
params.v            = 1;
params.beta         = 1;
params.tRate        = 0.01;
params.rhoSpines    = 1;
params.permeability = 1000;
params.phi          = 0.95;
params.maxRuns      = 10;
params.length       = 500;

% CaMKIIa-inspired parameters
params.D_m          = 0.001;
params.halflife_m   = 20*3600;
params.nucleotides  = 4300;

params.D_p          = 0.01; 
params.halflife_p   = 8*24*3600;
params.eta_p        = 32336;
params.aminoAcids   = 478;

%% distributions and costs
    
% iterate over strategies of interest
for stratInd = 1:2

    % apply either somatic or dendritic mRNA
    %----------------------------------------------------------------------    

    if stratInd == 1     % dendritic mRNA
        params.transp_m = 0.1;
        params.somaRet  = 0;
    elseif stratInd == 2 % somatic mRNA
        params.transp_m = 0;
        params.somaRet  = 1;
    end
    
    % fit ensemble mRNA diffusion coefficient 
    %----------------------------------------------------------------------    
    tmp                 = [];
    % MUST BE: vPlus = vMinus
    tmp.vPlus           = params.v;
    tmp.vMinus          = params.v;
    % MUST BE: bPlus = bMinus
    tmp.bPlus           = params.beta;
    tmp.bMinus          = params.beta;
    tmp.D               = params.D_m;
    tmp.transpProb      = params.transp_m;
    tmp.alpha           = tmp.transpProb/(1 - tmp.transpProb)*(1/(1/tmp.bPlus + 1/tmp.bMinus));
    tmp.halflife        = params.halflife_m;
    tmp.lambda          = log(2)/tmp.halflife;
    tmp.maxRuns         = params.maxRuns;
    % value at soma always set to 100
    tmp.atSoma          = 100;
    % for diffusion constant fit always set dendrite length to 100
    tmp.length          = 100;
    % solve 3-state model and fit 1-state model
    out                 = get_ensembleDiffConstFit(tmp);
    params.DEns_m       = out.D_1State;
    clear out tmp
    
    % no active protein transport
    params.transp_p     = 0;
    params.DEns_p       = params.D_p;

    % get mRNA and protein distributions
    %----------------------------------------------------------------------
    solution            = get_distribution(params);

    % Derive mRNA and protein distribution features
    %----------------------------------------------------------------------
    % initialise structures for mRNA and protein
    mRNA                = [];
    protein             = [];
    % Set the points at which the distributions should be evaluated. 
    xval                = linspace(0, params.length, 2*params.length);
    % Derive the mRNA and protein distributions. Round to 10 digits to 
    % remove errors below mass conservation accuracy that disturb the plots.
    mRNA.dendr           = round(deval(solution, xval, 2), 10);
    protein.dendr        = round(deval(solution, xval, 4), 10);
    % integrate dendritic mRNA and dendritic protein
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
    clear solution

    % get metabolic cost
    %----------------------------------------------------------------------
    cost                 = get_cost(params, mRNA, protein, xval);
    clear mRNA protein xval

    % store metabolic cost
    %----------------------------------------------------------------------    

    if stratInd == 1     % dendritic mRNA
        costDen = cost;
    elseif stratInd == 2 % somatic mRNA
        costSom = cost;
    end
    clear cost
end
clear stratInd

%% plot protein distributions and costs

plot_CaMK2aCost(costDen, costSom)

