function plot_somToDendrCostPerParam(data, varList, varsToShow)

% initialise figure
figure('FileName', [char(datetime('today', 'Format', 'yyyy_MM_dd')), '_somToDendrCostPerParam'], ...
       'Name'    , 'somToDendrCostPerParam', ...
       'Units'   , 'centimeter', ...
       'Position', [0, 0, 14, 10]);
% load colors
load('3Reds4Blues4Greens3Greys3Violets.mat');
% label of panels from top left to bottom right:
panelLabel         = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

% names of variables
varNameList           = varList.('Variable');
% sort data table along variables
data               = sortrows(data ,varNameList);

for index = 1:numel(varsToShow)
    tmpDat         = data;
    % variable name
    varName        = varNameList{varsToShow(index)};    
    % find variable values in the data table    
    varVals        = unique(tmpDat.(varName));    
    % number of values per parameter
    nVals          = numel(varVals);
    % check if exactly three values are present, otherwise give an error
    if nVals ~= 3
        error(['This plot requires exactly 3 values per parameter, but ', num2str(nVals), 'are present in the provided data.'])
    end
    % number of entries
    nRows          = size(tmpDat, 1);    
    % sort the data table along the values of the current variable of 
    % interest, then along all other parameters. Reshaping into "nVals"
    % columns then gives a similar ordering of parameter combinations
    % within each column
    sortOrder      = 1:size(varList, 1);
    sortOrder      = [varsToShow(index), sortOrder(sortOrder ~= varsToShow(index))];
    tmpDat         = sortrows(tmpDat, varNameList(sortOrder), 'ascend');
    clear sortOrder
    % get the logarithmic cost ratio of "somatic mRNA" vs 
    % "dendritic mRNA"
    somVsden       = tmpDat.('totalCost - som') ./ tmpDat.('totalCost - den');
    % reshape into matrix with entries for each parameter value in a single
    % column (works because table was sorted). The ordering of parameter 
    % combinations within each column is the same.
    somVsden       = reshape(somVsden, nRows/nVals, nVals);
    % Normalize with middle row, corresponding to the "middle" value of
    % the current parameter of interest. Note that this only works with 3
    % parameter values.
    somVsden(:, 1) = somVsden(:, 1)./somVsden(:, 2);
    somVsden(:, 3) = somVsden(:, 3)./somVsden(:, 2);
    somVsden(:, 2) = somVsden(:, 2)./somVsden(:, 2);
    % Finally, remove the baseline and transform to percent
    somVsden       = (somVsden - 1).*100;
    
    %% plot the results
    
    % plot three panels in the first and further panels in the second row
    if index < 3
        pos = [(index    )*4 + 1, 6, 2.5, 3];
    else
        pos = [(index - 3)*4 + 1, 1, 2.5, 3];
    end
    % major axis
    ax = axes('Units'               , 'centimeter', ...
              'Position'            , pos, ...
              'Color'               , 'none', ...
              'Box'                 , 'off', ...
              'FontSize'            , 10, ...
              'FontName'            , 'Arial', ...
              'LineWidth'           , 0.5, ...
              'TickDir'             , 'out', ...
              'TickLength'          , [0.02, 0.02], ...
              'TickLabelInterpreter', 'tex', ...
              'YAxisLocation'       , 'right', ...
              'XLim'                , [0.4, nVals + 0.6]); hold on;
    set(get(ax, 'XAxis'), 'Color', 'none')
    % panel label
    text('String'             , panelLabel{index}, ...
         'Interpreter'        , 'tex', ...
         'FontSize'           , 12, ...
         'FontName'           , 'Arial', ...
         'Units'              , 'centimeter', ...
         'Position'           , [-0.1, ax.Position(4)], ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment'  , 'bottom'); hold on;    
    % set the colors (reds for mRNA, colors.Blues for protein)
    if strcmp(varName, 'D_m') || strcmp(varName, 'halflife_m') || strcmp(varName, 'nucleotides:aminoAcids')
        shortColorList = [colors.LightRed; colors.Red; colors.DarkRed];
    elseif strcmp(varName, 'D_p') || strcmp(varName, 'halflife_p') || strcmp(varName, 'aminoAcids') || strcmp(varName, 'eta_p')
        shortColorList = [colors.LightBlue; colors.Blue; colors.DarkBlue];
    end
    % bar plot
    for bInd = 1:nVals
        yDat = median(somVsden(:, bInd));
        b = bar(bInd, yDat, ...
                'CData'    , shortColorList(bInd, :), ...
                'FaceColor', 'flat', ...
                'FaceAlpha', yDat > 0, ...
                'EdgeColor', 'flat', ...
                'LineWidth', 1); hold on;
        set(get(b, 'BaseLine'), 'Color', 'none')
    end
    % x-axis ticks
    for j = 1:nVals
        % correct units (hours for mRNA halflife, days for protein halflife,
        % log10 for diffusion constants)
        if strcmp(varName, 'halflife_m')
            tmpStr = num2str(varVals(j)/3600);
        elseif strcmp(varName, 'halflife_p')
            tmpStr = num2str(varVals(j)/(3600*24));
        elseif strcmp(varName, 'D_m')
            tmpStr = num2str(log10(varVals(j)));
        elseif strcmp(varName, 'D_p')
            tmpStr = num2str(log10(varVals(j)));
        elseif strcmp(varName, 'nucleotides:aminoAcids')
            tmpStr = [num2str(varVals(j) - 3), '/3'];
        else
            tmpStr = num2str(varVals(j));
        end
        text('String'             , tmpStr, ...
             'Interpreter'        , 'tex', ...
             'FontSize'           , 10, ...
             'FontName'           , 'Arial', ...
             'Units'              , 'normalized', ...
             'Position'           , [(2*j-1)/(2*nVals), 0], ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment'  , 'top'); hold on;
    end
    % x-axis label
    text('String'             , [varList.('Name for plots'){strcmp(varNameList, varName)}, ...
                                 varList.('Units for plots'){strcmp(varNameList, varName)}], ...
         'Interpreter'        , 'tex', ...
         'FontSize'           , 10, ...
         'FontName'           , 'Arial', ...
         'Units'              , 'normalized', ...
         'Position'           , [0.5, -0.2], ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment'  , 'middle'); hold on;
    % y-axis label
    if ismember(index, [2, 5])
        text('String'             , {'Somatic/dendritic cost', 'rel. to middle value [%]'}, ...
             'Interpreter'        , 'tex', ...
             'FontSize'           , 10, ...
             'FontName'           , 'Arial', ...
             'Units'              , 'normalized', ...
             'Position'           , [1.4, 0.5], ...
             'Rotation'           , -90, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment'  , 'bottom'); hold on;
    end
    % adjust y-axis limits
    yLim = get(gca, 'YLim');
    yLim = sign(yLim).*yLim;
    yLim = max(yLim);
    yLim = yLim.*[-1, 1];
    set(gca, 'YLim', yLim);
end

%% plot explanatory panel on the top-left figure edge

% axis
ax = axes('Units'               , 'centimeter', ...
          'Position'            , [2, 6, 1.5, 3], ...
          'Color'               , 'none', ...
          'Box'                 , 'off', ...
          'FontSize'            , 10, ...
          'FontName'            , 'Arial', ...
          'LineWidth'           , 0.5, ...
          'TickDir'             , 'out', ...
          'TickLength'          , [0.02, 0.02], ...
          'TickLabelInterpreter', 'tex', ...
          'XLim'                , [0, 1], ...
          'XTick'               , [], ...
          'YLim'                , [0, 1], ...
          'YTick'               , []); hold on;
set(get(ax, 'XAxis'), 'Color', 'none')
set(get(ax, 'YAxis'), 'Color', 'none')
% plot one upward-pointing and one downward-pointing arrow
plot([1, 1], [0, 1], ...
     'Color', 'black', ...
     'LineWidth', 0.5, ...
     'LineStyle', '-'); hold on;
plot(1, 0, ...
     'Color', 'black', ...
     'Marker', 'v', ...
     'MarkerSize', 3, ...
     'MarkerFaceColor', 'black'); hold on;
plot(1, 0.5, ...
     'Color', 'white', ...
     'Marker', 'o', ...
     'MarkerSize', 3, ...
     'MarkerFaceColor', 'white'); hold on;
plot(1, 1, ...
     'Color', 'black', ...
     'Marker', '^', ...
     'MarkerSize', 3, ...
     'MarkerFaceColor', 'black'); hold on;
% label upward-pointing arrow with "Dendritic mRNA gets more efficient"
text('String'             , {'Dendritic mRNA', 'preference'}, ...
     'Interpreter'        , 'tex', ...
     'FontSize'           , 10, ...
     'FontName'           , 'Arial', ...
     'Units'              , 'data', ...
     'Position'           , [0.9, 1], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment'  , 'top'); hold on;   
% label downward-pointing arrow with "Somatic mRNA gets more efficient"
text('String'             , {'Somatic mRNA', 'preference'}, ...
     'Interpreter'        , 'tex', ...
     'FontSize'           , 10, ...
     'FontName'           , 'Arial', ...
     'Units'              , 'data', ...
     'Position'           , [0.9, 0], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment'  , 'bottom'); hold on;  
