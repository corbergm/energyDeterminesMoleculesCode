function [sources, sourceDatSom, sourceDatNeur, signif] = get_proteinAbundanceVsNeuriteEnrichment()

%% Load data from Zappulo et al., 2017  

% read 'A) Differential expression' sheet
[~, ~, zappMRNAProt]               = xlsread('Zappulo_2017.xlsx', '2A) Differential expression');
% Load the gene names                                         in column "C"
%      the gene biotype                                       in column "D"
%      log2 protein in neurites                               in column "L"
%      log2 protein in somata                                 in column "M"
%      log2 RNA neurites:soma (RNA-seq)                       in column "T"

zappDat                            = cell2table(zappMRNAProt(4:34261, [3, 4, 12, 13, 20]));
zappDat.Properties.VariableNames   = {'gene name - Zappulo et al.', ...
                                      'gene biotype - Zappulo et al.', ...
                                      'protein quant. in neurites, log2 - Zappulo et al.', ...
                                      'protein quant. in somata, log2 - Zappulo et al.', ...
                                      'mRNA neurite:soma ratio, log2 - Zappulo et al.', ...
                                      };
clear zappMRNAProt

% get protein abundances for somata- and neurite enriched genes
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
% derive total protein abundance in neurites and soma
tmpSom  = zappDat.('protein quant. in somata, log2 - Zappulo et al.');
tmpNeur = zappDat.('protein quant. in neurites, log2 - Zappulo et al.');
nProt   = zeros(size(tmpSom));
for i = 1:size(tmpSom, 1)
    if ischar(tmpSom{i}) || ischar(tmpNeur{i})
        nProt(i) = NaN;
    else
        nProt(i) = 2.^tmpSom{i} + 2.^tmpNeur{i};
    end
end
clear i tmpNeur tmpSom
% transform back to log2 scale
nProt     = log2(nProt);
% keep only entries with valid enrichment scores and abundances
isValid   = ~isnan(neurToSom) & ~isnan(nProt);
neurToSom = neurToSom(isValid);
nProt     = nProt(isValid);
clear isValid
% get abundances for somata- and neurite-enriched entries
nProtSom  = nProt(neurToSom < -1);
nProtNeur = nProt(neurToSom >  1);
clear neurToSom nMRNA

% Wilcoxon rank sum test
% -------------------------------------------------------------------------

% compare protein abundances of somata-enriched and neurite-enriched entries
% via Wilcoxon rank sum test, check for significance
stars    = 0;
for alpha = [0.05, 0.01, 0.001, 0.0001]
    [~, h , ~] = ranksum(nProtSom, nProtNeur, 'alpha', alpha);
    if h == 1
        stars  = stars + 1;
    end
end
clear alpha h

%% Set data source and data

sources       = {'Zappulo et al., 2017'};
signif        = {stars};
sourceDatSom  = {nProtSom};
sourceDatNeur = {nProtNeur};
