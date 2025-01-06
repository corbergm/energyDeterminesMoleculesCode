function [sources, sourceDatSom, sourceDatNeur, signif] = get_mRNAAbundanceVsNeuriteEnrichment()

%% Load data from Zappulo et al., 2017  

% read 'A) Differential expression' sheet
[~, ~, zappMRNAProt]               = xlsread('Zappulo_2017.xlsx', '2A) Differential expression');
% Load the gene names                                         in column "C"
%      the gene biotype                                       in column "D"
%      log2 RNA in neurites (RNA-seq)                         in column "R"
%      log2 RNA in soma (RNA-seq)                             in column "S"
%      log2 RNA neurites:soma (RNA-seq)                       in column "T"
zappDat                            = cell2table(zappMRNAProt(4:34261, [3, 4, 18, 19, 20]));
zappDat.Properties.VariableNames   = {'gene name - Zappulo et al.', ...
                                      'gene biotype - Zappulo et al.', ...
                                      'mRNA quant. in neurites, log2 - Zappulo et al.', ...
                                      'mRNA quant. in somata, log2 - Zappulo et al.', ...
                                      'mRNA neurite:soma ratio, log2 - Zappulo et al.', ...
                                      };
clear zappMRNAProt

% get mRNA abundances for somata- and neurite enriched genes
% -------------------------------------------------------------------------

% derive mRNA neurite:soma enrichment
tmp       = zappDat.('mRNA neurite:soma ratio, log2 - Zappulo et al.');
neurToSom = zeros(size(tmp));
for i = 1:size(tmp, 1)
    if ischar(tmp{i})
        neurToSom(i) = NaN;
    else
        neurToSom(i) = tmp{i};
    end
end
clear i tmp
% derive total mRNA abundance in neurites and soma
nMRNA     = 2.^zappDat.('mRNA quant. in somata, log2 - Zappulo et al.') ...
          + 2.^zappDat.('mRNA quant. in neurites, log2 - Zappulo et al.');
% transform back to log2 scale
nMRNA     = log2(nMRNA);
% keep only entries with valid enrichment scores and abundances
isValid   = ~isnan(neurToSom) & ~isnan(nMRNA);
neurToSom = neurToSom(isValid);
nMRNA     = nMRNA(isValid);
clear isValid
% get abundances for somata- and neurite-enriched entries
nMRNASom  = nMRNA(neurToSom < -1);
nMRNANeur = nMRNA(neurToSom >  1);
clear neurToSom nMRNA

% Wilcoxon rank sum test
% -------------------------------------------------------------------------

% compare mRNA abundances of somata-enriched and neurite-enriched entries
% via Wilcoxon rank sum test, check for significance
stars    = 0;
for alpha = [0.05, 0.01, 0.001, 0.0001]
    [~, h , ~] = ranksum(nMRNASom, nMRNANeur, 'alpha', alpha);
    if h == 1
        stars  = stars + 1;
    end
end
clear alpha h

%% Set data source and data

sources       = {'Zappulo et al., 2017'};
signif        = {stars};
sourceDatSom  = {nMRNASom};
sourceDatNeur = {nMRNANeur};
