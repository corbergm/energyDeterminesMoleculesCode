

% author: CORNELIUS BERGMANN
% 
% last modified 03.02.2024

% This script calculates mRNA and protein distributions as well as the 
% associated costs for an exemplary set of parameters

% suppress solver detail output
warning('off', 'MATLAB:bvp5c:RelTolNotMet');

%% parameters

params              = [];

% mRNA parameters
params.D_m          = 0.001;     % Diffusion constant of non-transported dendritic mRNAs [micron^2/s]
params.halflife_m   = 20*3600;   % mRNA half-life [s]
params.nucleotides  = 4300;      % Number of nucleotides per transcript that have to be transcribed [a.u.]
params.somaRet      = 0.25;      % Share of mRNAs retained in the soma (from 0 to 1) [a.u.]
params.transp_m     = 0.1;       % Share of dendritic mRNAs that is transported at a given timepoint [a.u.]

% protein parameters
params.D_p          = 0.01;      % Diffusion constant of non-transported dendritic proteins [micron^2/s]
params.halflife_p   = 8*24*3600; % Protein half-life [s]
params.eta_p        = 1000;      % Maximal number of proteins in any spine [a.u.]
params.aminoAcids   = 478;       % Number of amino acids per protein [a.u.]
params.transp_p     = 0;         % Share of dendritic proteins that is transported at a given timepoint [a.u.]

% further parameters
params.v            = 1;         % default: 1     Instantaneous transport velocity [microns/s]
params.beta         = 1;         % default: 1     Switching rate out of transport state [1/s]
params.tRate        = 0.01;      % default: 0.01  Translation rate [protein/mRNA s]
params.rhoSpines    = 1;         % default: 1     Spine density along the dendrite [1/micron]
params.permeability = 1000;      % default: 1000  Effective ratio of proteins entering vs. leaving spines [a.u.]
params.phi          = 0.95;      % default: 0.95  Minimal filling rate of each spine (relative to the maximal filling eta_p) [a.u.]
params.maxRuns      = 10;        % default: 10    Maximal number of calls to the solver in "get_distribution" (not calls within "bvp5c")
params.length       = 500;       %                Dendrite length [microns]

%% get ensemble mRNA diffusion constant 
    
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

%% get ensemble protein diffusion constant 
    
tmp            = [];
% MUST BE: vPlus = vMinus
tmp.vPlus      = params.v;
tmp.vMinus     = params.v;
% MUST BE: bPlus = bMinus
tmp.bPlus      = params.beta;
tmp.bMinus     = params.beta;
tmp.D          = params.D_p;
tmp.transpProb = params.transp_p;
tmp.alpha      = tmp.transpProb/(1 - tmp.transpProb)*(1/(1/tmp.bPlus + 1/tmp.bMinus));
tmp.halflife   = params.halflife_p;
tmp.lambda     = log(2)/tmp.halflife;
tmp.maxRuns    = params.maxRuns;
% value at soma always set to 100
tmp.atSoma     = 100;
% for diffusion constant fit always set dendrite length to 100
tmp.length     = 100;
% solve 3-state model and fit 1-state model
out            = get_ensembleDiffConstFit(tmp);
params.DEns_p  = out.D_1State;
clear out tmp

%% get mRNA and protein distributions

solution = get_distribution(params);

%% Derive mRNA and protein distribution features

% initialise structures for mRNA and protein
mRNA                = [];
protein             = [];
% Set the points at which the distributions should be evaluated. 
xval                = linspace(0, params.length, 2*params.length);
% Derive the dendritic mRNA and protein distributions. Round to 10 digits
% (which is far below solver accuracy)
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
% derive absolute somatic mRNAs from the protein influx
mRNA.somaTotal      = protein.influx/params.tRate;
clear solution

%% get cost

cost = get_cost(params, mRNA, protein, xval);

%% plot mRNA and protein distributions and the associated costs

% load colors
load('3Reds4Blues4Greens3Greys3Violets.mat')

figure('WindowState', 'maximized')
% plot the mRNA distribution
subplot(1, 3, 1)
plot(xval, mRNA.dendr, 'Color', colors.Red)
xlabel('Dist. from soma [\mum]')
ylabel('mRNA')
set(gca, 'XLim'   , [0, params.length])
set(gca, 'Box'    , 'off')
set(gca, 'TickDir', 'out')
% plot the protein distribution
subplot(1, 3, 2)
plot(xval, protein.dendr + protein.spine, 'Color', colors.Blue)
xlabel('Dist. from soma [\mum]')
ylabel('Protein')
set(gca, 'XLim'   , [0, params.length])
set(gca, 'Box'    , 'off')
set(gca, 'TickDir', 'out')
% plot the associated transcription costs (cover mRNA synthesis and
% degradation), translation costs (cover protein synthesis and degradation)
% and transport costs (cover mRNA and potentially protein transport)
subplot(1, 3, 3)
bar(1:3, log10([cost.transcr, cost.transl, cost.transp]), 'FaceColor', colors.Grey)
ylabel('log_{10}cost [ATP/s]')
xticks(1:3)
xticklabels({'Transcr.', 'Transl.', 'Transp.'})
set(gca, 'Box'    , 'off')
set(gca, 'TickDir', 'out')
