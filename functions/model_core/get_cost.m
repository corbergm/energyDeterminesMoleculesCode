function cost = get_cost(params, mRNA, protein, xval)

% author: CORNELIUS BERGMANN
%
% modified 22.01.2023
%
% Input: 
% params   - parameter structure with fields
%            "v", "transp_m", "transp_p", "aminoAcids", "tRate", 
%            "nucleotides", "halflife_m", "length"
% mRNA     - mRNA data structure with fields
%            "dendr", "dendrTotal", "somaTotal"
% protein  - protein data structure with field
%            "dendr", "dendrTotal"
% xval     - x-values corresponding to "mRNA.dendr" and "protein.dendr"
%
% Output: 
% cost     - cost structure with fields
%            "transp", "transl", "transcr", "total", "totalPerSeg", "totalInSoma"
%

%% Output structure for metabolic cost

cost         = [];

%% Transport Cost

% Motor proteins use 125 ATP per second and the ratio of actively 
% transported mRNA (protein) is "transp_m" ("transp_p"), leading to transport 
% cost of
%
% C(transp)  = C(motor) v ( transp_m integral m + transp_p integral p )
%
cost.transp  = 125*params.v*(params.transp_m*mRNA.dendrTotal + params.transp_p*protein.dendrTotal);

%% Translation Cost

% Using 5 ATP per amino acid elongation, we get translation cost of
%
% C(transl)  = C(amino acid) n(amino acids) tau ( mSom + integral m )
%
cost.transl  = 5*params.aminoAcids*params.tRate*(mRNA.somaTotal + mRNA.dendrTotal);

%% Transcription Cost

% Using 2.17 ATP per nucleotide elongation, we get transcription cost of
%
% C(transcr) = C(nucleotide) n(nucleotide) lambda_m ( mSom + integral m )
%
cost.transcr = 2.17*params.nucleotides*log(2)/params.halflife_m*(mRNA.somaTotal + mRNA.dendrTotal);

%% Total Cost

cost.total   = cost.transp + cost.transl + cost.transcr;

%% Total cost per dendrite segments and in the soma

% length of segments
segLen           = 10;
% number of segments
nSegments        = ceil(params.length/segLen);
% container for total cost per dendritic segment
cost.totalPerSeg = zeros(1, nSegments);
% iterate over segments
for segInd = 1:nSegments
    % find indices within current segment
    inSeg     = xval >= (segInd - 1)*segLen & xval < segInd*segLen;
    % integrate dendritic mRNA and protein within current segment
    mDendrSeg = trapz(xval(inSeg), mRNA.dendr(inSeg));
    pDendrSeg = trapz(xval(inSeg), protein.dendr(inSeg));
    % integrate spine protein within current segment
    pSpineSeg = (protein.dendr(inSeg).*params.permeability)./((protein.dendr(inSeg).*params.permeability)./(params.rhoSpines*params.eta_p) + 1);
    pSpineSeg = trapz(xval(inSeg), pSpineSeg); 
    % total protein (dendrite + spine) within current segment
    pTotalSeg = pDendrSeg + pSpineSeg;
    % calculate total cost per segment, including ...
    costInSeg = 0;
    % ...translation cost:
    %     4 x dendritic mRNA x translation rate x protein length
    costInSeg = costInSeg ...
              + 4*mDendrSeg*params.tRate*params.aminoAcids;
    % ... degradation cost:
    %     1 x total protein x protein decay rate x protein length
    costInSeg = costInSeg ...
              + 1*pTotalSeg*log(2)/params.halflife_p*params.aminoAcids;
    % ... transport cost:
    %           125 x v x [ transp_m x dendritic mRNA ...
    %                     + transp_p x dendritic protein ]
    costInSeg = costInSeg ...
              + 125*params.v*(params.transp_m*mDendrSeg + params.transp_p*pDendrSeg);
    % store cost within current segment
    cost.totalPerSeg(segInd) = costInSeg;
end
% clear costInSeg inSeg nSegments mDendrSeg pDendrSeg pSpineSeg segInd
% cost in the soma include ...
cost.totalInSoma = 0;
% ...transcription cost:
%        2.17 x total mRNA x mRNA decay rate x mRNA length
cost.totalInSoma = cost.totalInSoma ...
                 + 2.17*(mRNA.dendrTotal + mRNA.somaTotal)*log(2)/params.halflife_m*params.nucleotides;
% ...translation cost:
%         4 x somatic mRNA x translation rate x protein length
cost.totalInSoma = cost.totalInSoma ...
                 + 4*mRNA.somaTotal*params.tRate*params.aminoAcids;
