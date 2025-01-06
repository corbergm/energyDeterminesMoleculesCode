function out = get_ensembleDiffConstFit(params)

% Input: 
% params   - parameter structure containing the fields:
%            "vPlus", "vMinus", "D", "length", "bPlus", "bMinus", "alpha",
%            "transpProb", "halflife", "atSoma", "maxRuns"
%
% Output: 
% out      - output structure with the fields
%            "D_1State"           (fitted ensemble diffusion constant)
%            "err_3State"         (error of 3-state model soultion)
%            "nMeshPts_3State"    (# mesh points used to solve 3-state model)
%            "anteroToRetro"      (array with antero/retrograde ratio)
%            "anteroRetroToTotal" (array with (antero+retrograde)/total ratio)
%            "err_fitting"        (error of fitting 1-state model to 3-state model)
%

% output structure
out                       = [];

%% solve 3-state model

% 3-state steady state model with
%                           a = alpha,
%                           b = bPlus = bMinus,
%                           v = vPlus = vMinus.
%
% (EQ I)   0 = - v yP' - b yP + a y0
% (EQ II)  0 =   v yM' - b yM + a y0
% (EQ III) 0 = b yP + b yM + D y0'' - 2 a y0 - lambda y0
%
% (BC I)   0 = yP (L) - yM(L)
% (BC II)  0 = y0'(L)
% (BC III) 0 = yP (0) - yM(0)
% (BC IV)  0 = y0 (0) + yP(0) + yM(0) - atSoma

% get solution
solution                  = get_steadyStateSol3stateModel(params);
% x-values to fit analytical 1-state solution
xvalForFit                = 0:(params.length/500):params.length;
% add up y0, yP, and yM to obtain the total concentration
anterograde_3State        = deval(solution, xvalForFit, 1);
retrograde_3State         = deval(solution, xvalForFit, 2);
resting_3State            = deval(solution, xvalForFit, 4);
total_3State              = anterograde_3State ...
                          + retrograde_3State ...
                          + resting_3State;
% store antero-to-retrograde and (antero-plus-retrograde)-to-total ratios
out.anteroToRetro         = anterograde_3State./retrograde_3State;
out.anteroRetroToTotal    = (anterograde_3State + retrograde_3State)./total_3State;
% store solver error and number of mesh points used
out.err_3State            = solution.maxErr;
out.nMeshPts_3State       = solution.nMeshPts;
clear anterograde_3State resting_3State retrograde_3State solution

%% fit analytical solution of 1-state model

% 1-state steady state model:
%
% (EQ I)   0 = D y'' - lambda y
%
% (BC I)   0 = y'(L)
% (BC II)  0 = y(0) - atSoma

% set least-squares algorithm tolerance
options                   = optimset('TolFun', 1e-8, 'TolX', 1e-8);
% fit 1-state model with ensemble diffusion constant "D_1State" and somatic 
% concentration "atSoma_1State" against 3-state model
[optVals, ~, exitflag, ~] = fminsearch(@(y)sum((get_analyticalSteadyStateSol1StateModel(y(1), xvalForFit, params.lambda, params.length, y(2)) - total_3State).^2), [0.1, params.atSoma], options);
out.D_1State              = optVals(1);
atSoma_1State             = optVals(2);
% check if fitting did work out
if exitflag ~= 1
    warning('Ensemble diffusion coefficient not fitted correctly')
end
clear exitflag options optVals

%% store results

% calculate and store L1-error of the fitted curve
total_1State              = get_analyticalSteadyStateSol1StateModel(out.D_1State, xvalForFit, params.lambda, params.length, atSoma_1State);
err_fitting               = trapz(xvalForFit, abs(total_1State - total_3State));
% to compare among multiple lengths, normalize with area under the
% curve of 3-state model
out.err_fitting           = err_fitting/trapz(xvalForFit, total_3State);
