

% author: CORNELIUS BERGMANN
% 
% last modified 12.07.2022

% This script calculates mRNA and protein distributions for CaMKIIa 
% and plots the predicted distributions as well as experimental data from
% Fonkeu et al., 2019 (Neuron).

% suppress solver detail output
warning('off', 'MATLAB:bvp5c:RelTolNotMet');

%% CaMK2a parameters for mRNA and protein

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
params.somaRet      = 0.25; 
params.transp_m     = 0.1;

params.D_p          = 0.01;      
params.halflife_p   = 8*24*3600; 
params.eta_p        = 32336;     
params.aminoAcids   = 478;
params.transp_p     = 0;

%% get CaMK2a mRNA and protein distributions
    
% fit ensemble mRNA diffusion coefficient 
%----------------------------------------------------------------------    
tmp            = [];
% MUST BE: vPlus = vMinus
tmp.vPlus      = params.v;
tmp.vMinus     = params.v;
% MUST BE: bPlus = bMinus
tmp.bPlus      = params.beta;
tmp.bMinus     = params.beta;
tmp.D          = params.D_m;
tmp.transpProb = params.transp_m;
tmp.alpha      = tmp.transpProb/(1 - tmp.transpProb)*(1/(1/tmp.bPlus + 1/tmp.bMinus));
tmp.halflife   = params.halflife_m;
tmp.lambda     = log(2)/tmp.halflife;
tmp.maxRuns    = params.maxRuns;
% value at soma always set to 100
tmp.atSoma     = 100;
% for diffusion constant fit always set dendrite length to 100
tmp.length     = 100;
% solve 3-state model and fit 1-state model
out            = get_ensembleDiffConstFit(tmp);
params.DEns_m  = out.D_1State;
clear out tmp

% no active protein transport
%----------------------------------------------------------------------
params.DEns_p = params.D_p;

% get mRNA and protein distributions
%----------------------------------------------------------------------
solution = get_distribution(params);

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
clear solution

% mRNA and protein distributions to be shown
mRNASim.x = xval;
mRNASim.y = mRNA.dendr;
protSim.x = xval;
protSim.y = protein.dendr + protein.spine;
clear xval mRNA protein

%% load protein and mRNA densities from Fonkeu et al., 2019

load('data\mRNADensity100Microns_Fonkeu_2019.mat')
mRNAExp.x = mRNADens(:, 1);
mRNAExp.y = mRNADens(:, 2);
clear mRNADens

load('data\proteinDensity100Microns_Fonkeu_2019.mat')
protExp.x = protDens(:, 1);
protExp.y = protDens(:, 2);
clear protDens

%% plot mRNA and protein distributions

plot_CaMK2aDistributions(params, mRNASim, protSim, mRNAExp, protExp)
