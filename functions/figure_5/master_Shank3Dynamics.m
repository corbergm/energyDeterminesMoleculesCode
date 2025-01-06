
% This script simulates the spatio-temporal spread of somatically
% photoactivated Shank3 proteins in a linear dendrite. Obtained results are
% stored under a defined date string "dateStr", for available reuslts see
% below. For plotting, one set of simulation results is loaded and plotted
% together with experimental data from Tsuriel et al., 2006.

%% parameters for Shank3/ProSAP2  

% minimal supply ratio
params              = [];
params.length       = 750; 
% ensemble mRNA diffusion coefficient, inludes active mRNA transport
params.DEns_m       = 0.1;    
% ensemble protein diffusion coefficient, no active protein transport
params.DEns_p       = 0.9;    
params.tRate        = 0.01;
params.halflife_m   = 12.68*3600;
params.halflife_p   = 20*24*3600; 
params.lambda_p     = log(2)/params.halflife_p; 
params.rhoSpines    = 1;
params.eta_p        = 50;
params.phi          = 0.95;      
params.somaRet      = 0.5; 

%% spine dynamics  

% rate of spine exit: 10 minutes
params.exitRate     = log(2)/1200; 
% rate of spine entering: 1000 x exit rate
params.uptakeRate   = 1000*params.exitRate;
% permeability = uptake / (exit + decay)
params.permeability = params.uptakeRate/(params.exitRate + params.lambda_p);

%% spatial and temporal grid step sizes

% Temporal grid with step 0.01s (10ms) and time limit 8h, unit is seconds:
DeltaT              = 0.01;
maxT                = 8*3600;
% Spatial grid on [0, L] with step 0.5microns, unit is microns:
DeltaX              = 0.5;

%% run simulation and save results

params.maxRuns      = 10;

% data string
dataStr             = '2023_05_04';
% explicit or implicit solution
exp_imp             = 'implicit';
% results are saved under ['\files\', dataStr, '_Shank3Dynamics_results_', exp_imp, '.mat']
run_Shank3Dynamics(exp_imp, dataStr, params, DeltaT, maxT, DeltaX)
