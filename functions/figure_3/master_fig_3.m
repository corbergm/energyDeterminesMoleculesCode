
%% This script will generate the data analysis in Figure 3
% Perform the statistical tests, and visualize the results
% ========================================================== Kanaan Mousaei
%

%% ***** Initiate script and make legend
close all
clear all
clc

fig3_labels = ['ABCDE'];
f = figure;
% f = openfig(".\Legend\legends.fig");
f.Units = 'centimeter';
f.Position = [10, 10, 13, 10];
pos = get(gca, 'Position');
pos(1) = 0.5; pos(2) = 6;
set(gca, 'Position', pos);

x_axis = [6, 10, 2, 6, 10];
y_axis = [6, 6, 1, 1, 1];

%% Run Figure 3 panel A

clearvars -except f x_axis y_axis fig3_labels

Tushev.data_location       = ".\data\Tushev_2018.xls";
Tushev.Sheet_name          = "PASSData";
% Load the data
Tushev.data                = read_data(Tushev.data_location,...
                                        Tushev.Sheet_name);
% Make the header titles for new list. Create filter properties for
% analyzing the reported data

Properties.Headers                  = [Tushev.data(1,:)];
Properties.filters.gene_feature     = '3pUTR';
Properties.filters.cell_type        = 'neuron-enriched';
Properties.filters.halflife_range   = [0, 24];
Properties.filters.rSquareMin       = 0.5;

% Filter the data using given filters and separate them into enrichment
% categories of somata or neuropil enriched
mrna                        = localization_tushev(Tushev.data,Properties);
% Extract only the half-life of each enrichment category
Halflife_localized.Somata               = double(mrna.Somata(2:end,14));
Halflife_localized.Neuropil             = double(mrna.Neuropil(2:end,14));
% Figure properties to plot and visualization
Fproperties.ccode           = 'm';
Fproperties.panel           = fig3_labels(1);
Fproperties.x_label         = ["Tushev et al., 2018"];
Fproperties.y_label         = "mRNA half-life [h]";

[Halflife_localized.Somata_processed,...
    Halflife_localized.Neuropil_processed,...
    condition.log,...
    condition.btstrp]       = operation(Halflife_localized.Somata,...
                                        Halflife_localized.Neuropil,...
                                        3);
% Perform statistical test on two localization category
[stat.p,stat.h]       =   ranksum(Halflife_localized.Somata,...
                                    Halflife_localized.Neuropil);
% boxplot the analysis with given properties
boxplot_full(Halflife_localized.Somata_processed,...
            Halflife_localized.Neuropil_processed,...
            Fproperties.ccode,...
            Fproperties.panel,...
            Fproperties.y_label,...
            Fproperties.x_label,stat.p,x_axis(1),y_axis(1))

%% Run Figure 3 panel B

clearvars -except f x_axis y_axis fig3_labels

Zappulo.data_location       = ".\data\Zappulo_2017.xlsx";
Zappulo.Sheet_name          = "2A) Differential expression";
% load the data
Zappulo.data                = read_data(Zappulo.data_location,...
                                        Zappulo.Sheet_name);
% Separate data into enrichment categories by dividing factor
Enrichment          = localization_zappulo(Zappulo.data,1);
% extract only protein length of each category and process them
% Transcript - CDS
% UTR.Somata      = (double(Enrichment.Somata(2:end,8))) - (double(Enrichment.Somata(2:end,9)));
% UTR.Neuropil    = (double(Enrichment.Neuropil(2:end,8))) - (double(Enrichment.Neuropil(2:end,9)));
UTR.Somata      = (double(Enrichment.Somata(2:end,8)));
UTR.Neuropil    = (double(Enrichment.Neuropil(2:end,8)));

% UTR.Somata(UTR.Somata==0)=[];
% UTR.Neuropil(UTR.Neuropil==0)=[];
% Perform statistical test on two localization category
[stat.p,stat.h]       =   ranksum(UTR.Somata,...
                                    UTR.Neuropil);
% Figure properties to plot and visualization
Fproperties.ccode           = 'm';
Fproperties.panel           = fig3_labels(2);
Fproperties.x_label         = ["Zappulo et al., 2017"];
Fproperties.y_label         = "mRNA length";

[UTR.Somata_processed,...
    UTR.Neuropil_processed,...
    condition.log,...
    condition.btstrp]       = operation(UTR.Somata,...
                                        UTR.Neuropil,...
                                        3);
% boxplot the analysis with given properties
boxplot_full(UTR.Somata_processed,...
            UTR.Neuropil_processed,...
            Fproperties.ccode,...
            Fproperties.panel,...
            Fproperties.y_label,...
            Fproperties.x_label,stat.p,x_axis(2),y_axis(2))

%% Run Figure 3 panel C

clearvars -except f x_axis y_axis fig3_labels Zappulo

Helm.data_location          = ".\data\Helm_2021.xlsx";
Helm.Sheet_name             = "Sheet1";
Helm.data                   = read_data(Helm.data_location,...
                                        Helm.Sheet_name);

Gene_ref.data_location      = ".\data\protein_gene_name.xlsx";
Gene_ref.Sheet_name         = "Name";
Gene_ref.data               = read_data(Gene_ref.data_location,...
                                        Gene_ref.Sheet_name);

load(".\data\2022_07_01_mouse_rat_homologues_MGI.mat");
Gene_ref.Translation        = mouseRatHomol;
clear mouseRatHomol
% Match found gene names with Helm data and translation refference of
% mouse/rat
Helm_genes                  = gene_name_match(Gene_ref.data,Helm.data);
Helm_genes.Translation      = gene_translation_match(Helm_genes.Matched,Gene_ref.Translation);
% Create header titles for matching Helm data (with gene names) to Zappulo
Headers.full                = ["Zappulo et al. 2017",...
                                Zappulo.data(3,:),...
                                Helm_genes.Translation(1,:)];
Headers.interested          = [Helm_genes.Translation(1,5:6),...
                                Helm_genes.Translation(1,3:4),...
                                Helm_genes.Translation(1,8),...
                                Helm_genes.Translation(1,11),...
                                "Zappulo et al. 2017",...
                                Zappulo.data(3,3),...
                                Zappulo.data(3,20),...
                                Zappulo.data(3,24),...
                                "Translation",...
                                Helm_genes.Translation(1,25:26)];
% Search, and match data of Helm and Zappulo. Then separate them into two
% categories of Neuropil-enriched and Somata-enriched
protein                     = cross_match_helm_zappulo(Helm_genes.Translation,Zappulo.data,Headers);
% extract only copy numbers of each category and process them
Copy_number.Somata              = text_extract(protein.Somata_interested);
Copy_number.Neuropil            = text_extract(protein.Neuropil_interested);
% Perform statistical test on two localization category
[stat.p,stat.h]       =   ranksum(Copy_number.Somata,...
                                    Copy_number.Neuropil);
% Figure properties to plot and visualization
Fproperties.ccode           = 'p';
Fproperties.panel           = fig3_labels(3);
Fproperties.x_label         = ["Helm et al., 2021",...
                                "Zappulo et al., 2017"];
Fproperties.y_label         = "Proteins/spine";
% boxplot the analysis with given properties
boxplot_full(Copy_number.Somata,...
            Copy_number.Neuropil,...
            Fproperties.ccode,...
            Fproperties.panel,...
            Fproperties.y_label,...
            Fproperties.x_label,stat.p,x_axis(3),y_axis(3))

%% Run Figure 3 panel D

clearvars -except f x_axis y_axis fig3_labels Zappulo

Fornasiero.data_location    = ".\data\Fornasiero_2018.xlsx";
Fornasiero.Sheet_name       = "Data";
Fornasiero.data             = read_data(Fornasiero.data_location,...
                                        Fornasiero.Sheet_name);
% Make the header titles for new list that matched between two database.
% One with shorters headers only to use parameters of interest. And one
% with full reported values of both database
Headers.full            = ["Fornasiero et al. 2018",...
                            Fornasiero.data(1,:),...
                            "Zappulo et al. 2017",...
                            Zappulo.data(3,:)];

Headers.interested      = [Fornasiero.data(1,3) + " (Fornasiero et al. 2018)",...
                            Fornasiero.data(1,4) + " (Fornasiero et al. 2018)",...
                            Zappulo.data(3,3) + " (Zappulo et al. 2017)",...
                            Zappulo.data(3,20) + " (Zappulo et al. 2017)",...
                            Zappulo.data(3,24) + " (Zappulo et al. 2017)",...
                            Fornasiero.data(1,5:14),...
                            Fornasiero.data(1,16:19)];
% Cross-match data and categorize into neuropil-enriched and
% somata-enriched
Enrichment             = cross_match_Fornasiero_zappulo(Fornasiero.data,...
                                        Zappulo.data,...
                                        1,...
                                        Headers);
% extract only half-lives of each category and process them
Halflife_localized.Somata            = double(Enrichment.Somata_interested(2:end,6));
Halflife_localized.Neuropil          = double(Enrichment.Neuropil_interested(2:end,6));

Halflife_localized.Somata(isnan(Halflife_localized.Somata))         = [];
Halflife_localized.Neuropil(isnan(Halflife_localized.Neuropil))     = [];
% Perform statistical test on two localization category
[stat.p,stat.h]       =   ranksum(Halflife_localized.Somata,...
                                    Halflife_localized.Neuropil);
% Figure properties to plot and visualization
Fproperties.ccode           = 'p';
Fproperties.panel           = fig3_labels(4);
Fproperties.x_label         = ["Fornasiero et al., 2018",...
                                "Zappulo et al., 2017"];
Fproperties.y_label         = "Protein half-life [d]";

[Halflife_localized.Somata_processed,...
    Halflife_localized.Neuropil_processed,...
    condition.log,...
    condition.btstrp]       = operation(Halflife_localized.Somata,...
                                        Halflife_localized.Neuropil,...
                                        3);
% boxplot the analysis with given properties
boxplot_full(Halflife_localized.Somata_processed,...
            Halflife_localized.Neuropil_processed,...
            Fproperties.ccode,...
            Fproperties.panel,...
            Fproperties.y_label,...
            Fproperties.x_label,stat.p,x_axis(4),y_axis(4))

%% Run Figure 3 panel E

clearvars -except f x_axis y_axis fig3_labels Zappulo

% Separate data into enrichment categories by dividing factor
Enrichment          = localization_zappulo(Zappulo.data,1);
% extract only protein length of each category and process them
PLength.Somata      = (double(Enrichment.Somata(2:end,9)))./3;
PLength.Neuropil    = (double(Enrichment.Neuropil(2:end,9)))./3;
% Perform statistical test on two localization category
[stat.p,stat.h]       =   ranksum(PLength.Somata,...
                                    PLength.Neuropil);
% Figure properties to plot and visualization
Fproperties.ccode           = 'p';
Fproperties.panel           = fig3_labels(5);
Fproperties.x_label         = ["Zappulo et al., 2017"];
Fproperties.y_label         = "Protein length";

[PLength.Somata_processed,...
    PLength.Neuropil_processed,...
    condition.log,...
    condition.btstrp]       = operation(PLength.Somata,...
                                        PLength.Neuropil,...
                                        3);
% boxplot the analysis with given properties
boxplot_full(PLength.Somata_processed,...
            PLength.Neuropil_processed,...
            Fproperties.ccode,...
            Fproperties.panel,...
            Fproperties.y_label,...
            Fproperties.x_label,stat.p,x_axis(5),y_axis(5))
%% Functions
function [data] = read_data(location,sheet)
    % Read data from given excel sheet, and convert it to string output
    temp_data       = readmatrix(location,"sheet",sheet);
    data_size       = size(temp_data);
    data_opts       = spreadsheetImportOptions("NumVariables",...
                                                data_size(2),...
                                                "Sheet",...
                                                sheet);
    data_temp       = readmatrix(location,...
                                  data_opts,...
                                  "UseExcel",...
                                  false);

    data            = string(data_temp);
end

function boxplot_full(x1,x2,c,pl,y_info,x_info,pv, ax1, ax2)
    % load colors
    DarkRed      = [0.9137, 0.4157, 0.4392];
    DarkBlue     = [0.4667, 0.6706, 0.7412];
    % proteins: blue, mRNA: red
    if strcmp(c,'m')
        color        = DarkRed;     % DarkRed or DakBlue
    elseif strcmp(c,'p')
        color        = DarkBlue;     % DarkRed or DakBlue
    end

    if length(x1(1,:))<length(x1(:,1))
        q={x1,x2};
    elseif length(x1(1,:))>length(x1(:,1))
        q={x1',x2'};
    end

    yNameList    = x_info;
    yLabelString = y_info;
    panelLabel   = pl;
    
    % create axis
    ax                             = axes;
    ax.Units                       = 'centimeter';
    ax.Position                    = [ax1, ax2, 2, 3];
    % plot panel label
    T                              = text; hold on;
    T.String                       = panelLabel;
    T.Interpreter                  = 'tex';
    T.FontSize                     = 12;
    T.FontName                     = 'Arial';
    T.Units                        = 'centimeter';
    T.Position                     = [-1.3, ax.Position(4) + 0.2];
    T.HorizontalAlignment          = 'center';
    T.VerticalAlignment            = 'bottom';
    
    for i = 1:2
        % plot boxplot ------------------------> here's the data

        if i == 1

            b                          = boxchart(i.*ones(size(q{i}, 1), 1), q{i}); hold on;
            b.BoxFaceColor             = color;
            b.BoxWidth                 = 0.5;
            b.BoxFaceAlpha             = 0;
            b.WhiskerLineColor         = color;
            b.LineWidth                = 1;
            b.MarkerColor              = 'none';

        elseif i == 2
            b                          = boxchart(i.*ones(size(q{i}, 1), 1), q{i}); hold on;
            b.BoxFaceColor             = color;
            b.BoxWidth                 = 0.5;
            b.BoxFaceAlpha             = 1;
            b.WhiskerLineColor         = color;
            b.LineWidth                = 1;
            b.MarkerColor              = 'none';

            p                     = plot(2 + (b.BoxWidth/2).*[-1, 1], median(q{i}).*[1, 1]); hold on;
            p.Color               = 'white';
            p.LineWidth           = 1;


        end
            set(gcf, 'InvertHardCopy', 'off');
            set(gcf, 'Color', 'w');   
    end
    
    % adjust axis layout
    ax.Box                         = 'off';
    ax.Color                       = 'none';
    ax.FontSize                    = 10;
    ax.FontName                    = 'Arial';
    ax.LineWidth                   = 1;
    ax.TickDir                     = 'out';
    ax.TickLength                  = [0.02, 0.03];
    ax.TickLabelInterpreter        = 'tex';
    ax.XAxis.Color                 = 'none';
    ax.XLim                        = [0.4, 2.6];
    ax.YLabel.String               = yLabelString;
    ax.YLabel.FontSize             = 10;
    ax.YLabel.FontName             = 'Arial';
    ax.YLabel.Interpreter          = 'tex';
    

    % plot x-axis tick labels
    T                          = text; hold on;
    T.String                   = yNameList;
    T.Interpreter              = 'tex';
    T.FontSize                 = 10;
    T.FontName                 = 'Arial';
    T.Units                    = 'normalized';
    T.Position                 = [0.5, -0.05];
    T.HorizontalAlignment      = 'center';
    T.VerticalAlignment        = 'top';




    if pv >= 0.05
        str     =   "ns";
        nStars  = 0;
    elseif 0.05 > pv && pv >= 0.01
        str     =   "*";
        nStars  = 1;
    elseif 0.01 >= pv && pv >= 0.001
        str     =   "**";
        nStars  = 2;
    elseif 0.001 > pv
        str     =   "***";
        nStars  = 3;
    end 

mat                            = [20,  0,  0,  0;
                                   0, 20,  0,  0;
                                   0,  0, 20,  0;
                                   0, 21,  0,  0]./20;
axS                            = axes('units', 'centimeter', 'Position', ax.Position*mat + [0, 0, 0, 0.3]);

% plot significance bar ...
p                              = plot([1, 2], [0.25, 0.25]); hold on;
p.Color                        = 'black';
p.LineWidth                    = 0.75;

T                              = text; hold on;
if nStars > 0
    T.String                   = repmat('\ast', 1, nStars);
else
    T.String                   = 'ns';
end
T.Interpreter                  = 'tex';
T.FontSize                     = 10;
T.FontName                     = 'Arial';
T.Units                        = 'data';
T.Position                     = [1.5, 0.25];
T.HorizontalAlignment          = 'center';
T.VerticalAlignment            = 'bottom';

% adjust axis layout
axS.Box                        = 'off';
axS.Color                      = 'none';
axS.XAxis.Color                = 'none';
axS.XLim                       = ax.XLim;
axS.YAxis.Color                = 'none';
axS.YLim                       = [0, 1];

end

function [out1,out2,out3,out4]=operation(inpt1,inpt2,cnd)
    % Performs additional visual processing on the data and output flags
    % condition 1 ==> original
    % condition 2 ==> log10
    % condition 3 ==> 10k bootstrap over mean
    % condition 4 ==> condition 2 then condition 3
    
    if cnd==1
        out1=inpt1;
        out2=inpt2;
        out3=0;
        out4=0;
    
    elseif cnd == 2
        
        out1=log10(inpt1);
        out2=log10(inpt2);
        out3=1;
        out4=0;
    
    elseif cnd ==3
        
        out1=bootstrp(10000,@mean,inpt1);
        out2=bootstrp(10000,@mean,inpt2);
        out3=0;
        out4=1;
        
    
    elseif cnd == 4
    
        out1_temp=log10(inpt1);
        out2_temp=log10(inpt2);
    
        out1=bootstrp(10000,@mean,out1_temp);
        out2=bootstrp(10000,@mean,out2_temp);
        
        out3=1;
        out4=1;
    end
end

function [output] = localization_tushev(input,props)
    % Create separate list with given header titles for each category
    output.Neuropil         = [props.Headers];
    output.Somata           = [props.Headers];

    for m = 2:length(input)
        % Filter the data will given parameters
        if strncmp(props.filters.gene_feature,input(m,3),5) && ...
                strcmp(props.filters.cell_type,input(m,8)) && ...
                double(input(m,14)) > min(props.filters.halflife_range) && ...
                double(input(m,14)) < max(props.filters.halflife_range) %&& ...
%                 double(input(m,8)) > props.filters.rSquareMin

            % Localize them into neuropil or somata enriched based on their
            % provided label in the dataset
            if strcmp(input(m,11),'somata')
                output.Somata           = [output.Somata;input(m,:)];
            elseif strcmp(input(m,11),'neuropil')
                output.Neuropil         = [output.Neuropil;input(m,:)];
            end
        end
    end
end

function [output] = localization_zappulo(input,dividing_factor)
    % Creat two categories of Somata-enriched and Neuropil-enriched
    output.Neuropil     = input(3,:);
    output.Somata       = input(3,:);
    
    % Divide data into two groups and eliminate missing data
    for p = 4:length(input)
        if ~ismissing(double(input(p,20))) && ~ismissing(double(input(p,8))) && ~ismissing(double(input(p,9)))
            if double(input(p,20)) > +abs(dividing_factor)
                % Neuropil-enriched category bigger than (+dividing favtor)
                output.Neuropil         = [output.Neuropil;input(p,:)];
            elseif double(input(p,20)) < -abs(dividing_factor)
                % Somata-enriched category less than (-dividing favtor)
                output.Somata           = [output.Somata;input(p,:)];
            end
        end
    end
end

function [output] = gene_name_match(input1,input2)
    % Creat a table and add matched gene names and Helm data all together.
    % (preprocessing for next step)
    output.Matched         = ["Manual",...
                                    input1(1,:),...
                                    "Helm et al. 2021",...
                                    input2(1,:)];

    % Find matching proteins in both table and add them correspondingly 
    for p = 2:length(input2)
        clear position
        position            = find(contains(input1(:,1),input2(p,1)));
        if ~isempty(position)
            for id = 1:length(position)
                if input2(p,1) == input1(position(id),1)
                    output.Matched          = [output.Matched;...
                                                    "*******",...
                                                    input1(position(id),:),...
                                                    "*******",...
                                                    input2(p,:)];
                end
            end
        elseif isempty(position)
            output.Matched              = [output.Matched;...
                                            "*******",...
                                            "","","yes",...
                                            "*******",...
                                            input2(p,:)];
        end
    end
end


function output = gene_translation_match(input1,input2)
    % Get mouse vs rat gene name translation data and match it to Helm data
    output          = [input1(1,:),...
                        "*******",...
                        "gene name (mouse)",...
                        "gene name (rat)"];
    translation_ref = string(table2cell(input2));
    for p = 2:length(input1)
        if ~strcmp(input1(p,4),"yes")
            clear gene_name
            gene_name                   = split(input1(p,3),", ");
            match_flag                  = 0;
            for s = 1:length(gene_name)
                name                = gene_name(s);
                position            = find(contains(lower(translation_ref(:,2)),lower(name)));
                if match_flag == 0
                    if isempty(position)
                        if gene_name(length(gene_name)) == gene_name(s)
                            output          = [output;...
                                                input1(p,:),...
                                                "*******",...
                                                "",...
                                                ""];
                        end
                    elseif ~isempty(position)
                        for id = 1:length(position)
                            if translation_ref(position(id),2) == name && match_flag == 0
                                output          = [output;
                                                    input1(p,:),...
                                                    "*******",...
                                                    translation_ref(position(id),:)];
                                match_flag = 1;
                            end
                        end
                    end
                elseif match_flag == 1
                    continue
                end
            end
        elseif strcmp(input1(p,4),"yes")
            output          = [output;...
                                input1(p,:),...
                                "*******",...
                                "",...
                                ""];
        end
    end
end


function [output] = cross_match_helm_zappulo(input1,input2,titles)
    % Creating tables with given full and interested header titles to match
    % and categories them into their enrichment of neuropil or somata
    output.Matched_full         = [titles.full];
    output.Matched_interested   = [titles.interested];
    output.Neuropil_full        = [titles.full];
    output.Neuropil_interested  = [titles.interested];
    output.Somata_full          = [titles.full];
    output.Somata_interested    = [titles.interested];
    
    % Avoid the proteins and genes that have subunits
    % And search the rest for matching genes and add them to corresponding
    % table
    for p = 2:length(input1)
        if strcmp(input1(p,4),"yes") || ismissing(input1(p,3))
            clear temp_list
            temp_list.full                = ["*******",...
                                            input2(1,:),...
                                            input1(p,:)];
            temp_list.interested          = [input1(p,5:6),...
                                            input1(p,3:4),...
                                            input1(p,8),...
                                            input1(p,11),...
                                            "*******",...
                                            input2(1,3),...
                                            input2(1,20),...
                                            input2(1,24),...
                                            "*******",...
                                            input1(p,25:26)];

            output.Matched_full         = [output.Matched_full;...
                                            temp_list.full];
            output.Matched_interested   = [output.Matched_interested;...
                                            temp_list.interested];
            continue
        elseif ~strcmp(input1(p,4),"yes")
            clear gene_name
            gene_name               = split(input1(p,3),", ");
            match_flag              = 0;
            for s = 1:length(gene_name)
                name                = gene_name(s);
                position            = find(contains(lower(input2(:,3)),lower(name)));
                if match_flag == 0
                    if ~isempty(position)
                        for id = 1:length(position)
                            if input2(position(id),3) == name && match_flag == 0
                                clear temp_list
                                temp_list.full                = ["*******",...
                                                                input2(position(id),:),...
                                                                input1(p,:)];
                                temp_list.interested          = [input1(p,5:6),...
                                                                input1(p,3:4),...
                                                                input1(p,8),...
                                                                input1(p,11),...
                                                                "*******",...
                                                                input2(position(id),3),...
                                                                input2(position(id),20),...
                                                                input2(position(id),24),...
                                                                "*******",...
                                                                input1(p,25:26)];

                                output.Matched_full         = [output.Matched_full;...
                                                                temp_list.full];
                                output.Matched_interested   = [output.Matched_interested;...
                                                                temp_list.interested];
    
                                match_flag = 1;
                            end
                        end
                    end
                elseif match_flag == 1
                    continue
                end
            end
        end
    end

    % Localize data into somata or neuropil enriched category and avoid
    % subunits and check their rat vs mouse translation names match
    for localization = 2:length(output.Matched_full(:,1))
        if ~strcmp(output.Matched_full(localization,33),"yes") || strcmp(output.Matched_full(localization,54),output.Matched_full(localization,55))
            if ~strcmp(output.Matched_full(localization,55),"")
                if double(output.Matched_full(localization,21)) > 0
                    % neuropil-enriched category bigger than 0
                    output.Neuropil_full            = [output.Neuropil_full;...
                        output.Matched_full(localization,:)];
                    output.Neuropil_interested      = [output.Neuropil_interested;...
                        output.Matched_interested(localization,:)];
                elseif double(output.Matched_full(localization,21)) < 0
                    % Somata-enriched category less than 0
                    output.Somata_full              = [output.Somata_full;...
                        output.Matched_full(localization,:)];
                    output.Somata_interested        = [output.Somata_interested;...
                        output.Matched_interested(localization,:)];
                end
            end
        end
    end
end

function [output] = text_extract(input)
    % extract main measurement of copy number from string cell
    output = [];
    for cp = 2:length(input)
        clear temp_cp
        temp_cp = split(input(cp,6));
        output = [output;temp_cp(1)];
    end
    output = double(output);
end

function [output] = cross_match_Fornasiero_zappulo(input_data1,input_data2,dividing_factor,titles)
    % Run the algorithm that searches gene names in both database and find
    % the matching ones. Then create new table with matched genes and
    % separate them into localization (neuropil or somata enriched) based
    % on the mRNA enrichement scores.

    % Defining new categories and their header titles
    output.df                   = dividing_factor;
    output.Matched_full         = [titles.full];
    output.Matched_interested   = [titles.interested];
    output.Neuropil_full        = [titles.full];
    output.Neuropil_interested  = [titles.interested];
    output.Somata_full          = [titles.full];
    output.Somata_interested    = [titles.interested];

    % Search genes in both database and find matching ones, put them into
    % mathced list
    for row = 2:length(input_data1)
        clear gene_name
        gene_name               = split(input_data1(row,3));
        match_flag              = 0;
        for s = 1:length(gene_name)
            name                = gene_name(s);
            position            = find(contains(lower(input_data2(:,3)),lower(name)));
            if match_flag == 0
                if ~isempty(position)
                    for id = 1:length(position)
                        if input_data2(position(id),3) == name
                            clear temp_list
                            temp_list.full       = ["Fornasiero",...
                                                    input_data1(row,:),...
                                                    "Zappulo",...
                                                    input_data2(position(id),:)];
                            temp_list.interested = [input_data1(row,3),...
                                                    input_data1(row,4),...
                                                    input_data2(position(id),3),...
                                                    input_data2(position(id),20),...
                                                    input_data2(position(id),24),...
                                                    input_data1(row,5:14),...
                                                    input_data1(row,16:19)];

                            output.Matched_full         = [output.Matched_full;...
                                                            temp_list.full];
                            output.Matched_interested   = [output.Matched_interested;...
                                                            temp_list.interested];

                            match_flag = 1;
                        end
                    end
                end
            elseif match_flag == 1
                continue
            end
        end
    end

    % Separate genes based on their mRNA enrichment and localization
    % according to provided dividing factor
    for localization = 2:length(output.Matched_full(:,1))
        if double(output.Matched_full(localization,109)) > +(abs(dividing_factor))
            % Neuropil-enriched category bigger than (+ dividing factor)
            output.Neuropil_full            = [output.Neuropil_full;...
                output.Matched_full(localization,:)];
            output.Neuropil_interested      = [output.Neuropil_interested;...
                output.Matched_interested(localization,:)];
        elseif double(output.Matched_full(localization,109)) < -(abs(dividing_factor))
            % Somata-enriched category less than (- dividing factor)
            output.Somata_full              = [output.Somata_full;...
                output.Matched_full(localization,:)];
            output.Somata_interested        = [output.Somata_interested;...
                output.Matched_interested(localization,:)];
        end

    end
end