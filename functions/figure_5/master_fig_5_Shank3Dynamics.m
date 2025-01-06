
%% load simulation results

% with the explicit or the implicit finite difference scheme
% load('files\2023_05_04_Shank3Dynamics_results_explicit');
load('files\2023_05_04_Shank3Dynamics_results_implicit');

%% load data from Tsuriel et al., 2006

% DOI: 10.1371/journal.pbio.0040271
% Data were extracted manually from Figure 8D of (Tsuriel et al., 2006)
% using WebPlotDigitizer (https://automeris.io/WebPlotDigitizer).

% To allow normalization of each curve, we binned the data with a time
% step of 0.2h and averaged values falling in one bin. 
load('data\Shank3Dynamics_binned_Tsuriel_2006')

% store data in one matrix
datTsuriel       = zeros(4, 39);
datTsuriel(1, :) = sp25;
datTsuriel(2, :) = sp50;
datTsuriel(3, :) = sp75;
datTsuriel(4, :) = sp105;
% use the same markers as (Tsuriel et al., 2006)
spineMarker      = {'d', '^', 'o', 's'};
% use the same average spine locations as  (Tsuriel et al., 2006)
spineLocs        = [25, 50, 75, 105];
clear sp25 sp50 sp75 sp105

%% plot results

% times to plot Shank3 distributions 8h, 4h, 1h after photoactivation onset
timeIndsToPlot = [17, 9, 3];
plot_Shank3DynamicsAndTsuriel2006Data(params, x, pSom, pInit, pInitSpine, pTag, pTagSpine, t, timeIndsToPlot, saveTimes, spineLocs, spineMarker, measuredTimes, datTsuriel)

clear dateStr datTsuriel f measuredTimes panelLabel params pInit pInitSpine 
clear pSom pTag pTagSpine saveTimes spineLocs spineMarker t timeIndsToPlot x
