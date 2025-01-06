
%% load data and prepare data

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
    % remove cost per 50 microns, table concatenation for multiple length
    % does not allow it
    datTmp   = removevars(datTmp,{'totalCostPer10MicronSegment - som', 'totalCostPer10MicronSegment - den'});
    % Append data
    if indDate == 1
        data = datTmp;
    else
        data = [data; datTmp];
    end
end
clear dataStr datTmp indDate

%% plot total cost distribution vs protein count

% dendrite length of interest
lengthToShow = 500;
% y-limit of plots
yLim         = [-1, 10];

plot_globalCost(data(data.('length') == lengthToShow, :), yLim)
clear lengthToShow yLim
