function run_Shank3Dynamics(exp_imp, dateStr, params, DeltaT, maxT, DeltaX)

% This function simulates the spatio-temporal spread of somatically
% photoactivated Shank3 protein with parameters "params" on a spatial
% (temporal) grid "x" ("t") with step size "DeltaX" ("DeltaT")
%
% Input: 
% params   - parameter structure containing the fields:
%            "length", "DEns_m", "DEns_p", "tRate", "halflife_m", 
%            "halflife_p", "lambda_p", "rhoSpines", "eta_p", "phi", 
%            "somaRet", "exitRate", "uptakeRate", "permeability", "maxRuns"
% DeltaT   - temporal step size in seconds
% maxT     - time to simulate
% DeltaX   - spatial step size in microns
%
% Note: results are saved under ['\files\Shank3Dynamics\', dateStr, '_results.mat']

%% spatial and temporal grid  

% Temporal grid with step "DeltaT" and time limit "maxT", unit is seconds:
t                             = 0:DeltaT:maxT;      
% Spatial grid on [0, L] with step "DeltaX", unit is microns:
x                             = 0:DeltaX:params.length; 
% size of spatial grid
N                             = numel(x); 

%% steady state distribution for non-tagged proteins  

solution                      = get_distribution(params);
% steady state somatic concentration of non-tagged proteins
pSom                          = deval(solution, 0, 4);
% steady state distribution of non-tagged protein
pInit                         = deval(solution, x, 4);
pInitSpine                    = params.permeability.*pInit./(params.permeability/(params.rhoSpines*params.eta_p).*pInit + 1);
% store spine-to-dendrite ratio (division of sums is sufficient because
% "pInit" and "pInitSpine" are defined on the same x-grid, no integration
% needed here)
params.spineToDendriteRatio   = sum(pInitSpine)/sum(pInit);

% for later use, define the factor 
%
%                    pInitSpine
% supplyRatio = 1 - ------------
%                     rho eta  
%
supplyRatio                   = 1 - pInitSpine./(params.rhoSpines*params.eta_p);

% figure;plot(1 - supplyRatio)
% figure; plot(solution.x, solution.y(2,:));
figure; subplot(1, 2, 1); plot(pInit); subplot(1, 2, 2); plot(pInitSpine);

clear solution

%% spatio-temporal evolution of photoactivated proteins  

% time point indices where to save distributions (at 0:0.5:8 hours)
saveTimes       = find(ismember(t, 3600.*(0:0.5:8)));
% tagged, fluorescent proteins, measured at 100 time points
pTag            = zeros(numel(saveTimes), N);
pTagSpine       = zeros(numel(saveTimes), N);
% save first time step distribution (all zeros)
pTag(1, :)      = zeros(1, N);
pTagSpine(1, :) = zeros(1, N);

% increase saveInd
saveInd         = 2;
   
% use either an explicit or implicit finite difference method
if strcmp(exp_imp, 'explicit')

    % coefficient matrix for explicit method
    % ---------------------------------------------------------------------
    
    % matrix size: N x N

    %  | p(1, i+1)  |     | 1 0  0            | 0               |     | p(1, i)  |
    %  | p(2, i+1)  |     | b a1 b            |   c             |     | p(2, i)  |
    %  | p(3, i+1)  |     |   b a2 b          |     c           |     | p(3, i)  |
    %  |    ...     |     |     .  . .        |       .         |     |   ...    |
    %  |    ...     |     |       .  . .      |         .       |     |   ...    |
    %  |    ...     |     |         .  . .    |           .     |     |   ...    |
    %  |    ...     |     |           b  . b  |             c   |     |   ...    |
    %  | p(N, i+1)  |     |             0  aN |               c |     | p(N, i)  |
    %  |            |  =  |-------------------------------------|  x  |          |
    %  | pS(1, i+1) |     | d1                | c               |     | pS(1, i) |
    %  | pS(2, i+1) |     |   d2              |   c             |     | pS(2, i) |
    %  | pS(3, i+1) |     |     d3            |     c           |     | pS(3, i) |
    %  |    ...     |     |        .          |       .         |     |   ...    |
    %  |    ...     |     |          .        |         .       |     |   ...    |
    %  |    ...     |     |            .      |           .     |     |   ...    |
    %  |    ...     |     |              .    |             c   |     |   ...    |
    %  | pS(N, i+1) |     |                dN |               c |     | pS(N, i) |
    
    %           2 D_p dT
    % aK = 1 - ---------- - lambda dT - supplyRatio(xK) uptakeRate dT
    %             dX^2
    %
    %       D dT
    % b  = ------
    %       dX^2
    %  
    % c  = exitRate dT
    %
    % d  = supplyRatio(xK) uptakeRate dT
    %
    % e  = 1 - exitRate dT - lambda dT 
    
    % set up coefficient matrix for explicit method
    % ---------------------------------------------------------------------

    % left upper
    a     = 1 - 2*params.DEns_p*DeltaT/DeltaX^2 - params.lambda_p*DeltaT - supplyRatio.*params.uptakeRate.*DeltaT;
    b     = params.DEns_p*DeltaT/DeltaX^2;
    b     = b.*ones(1, N-2);
    LU    = diag([b, 0], -1) + diag([1, a(2:N)], 0) + diag([0, b], 1);
    % right upper
    c     = params.exitRate*DeltaT;
    c     = [0, c.*ones(1, N-1)];
    RU    = diag(c, 0);
    % left lower
    d     = supplyRatio.*params.uptakeRate.*DeltaT;
    LL    = diag(d, 0);
    % right lower
    e     = 1 - params.exitRate*DeltaT - params.lambda_p*DeltaT;
    e     = e.*ones(1, N);
    RL    = diag(e, 0);
    % consolidate in one matrix
    A_exp = sparse([LU, RU; LL, RL]);
    % check that all entries are within [-1, 1]
    if 0 < sum(abs(a) >= 1) + sum(abs(b) >= 1) + sum(abs(c) >= 1) + sum(abs(d) >= 1) + sum(abs(e) >= 1)
        warning('Coefficient matrix for explicit method contains entries bigger than 1.')
    end
    clear a b c d e LL LU RL RU

    % run explicit method
    % ---------------------------------------------------------------------

    % counter for progress of script
    done                          = 0;
    fprintf('\n %d percent done', done);

    in    = sparse(zeros(2*N, 1));
    in(1) = pSom;

    for tInd = 2:numel(t)
        % print progress 
        if round(tInd/numel(t)*100) > done
            done                  = round(tInd/numel(t)*100);
            fprintf('\n %d percent done', done);
        end
        % multiply with coefficient matrix
        out = A_exp*in;
        % re-set input array
        in  = out;
        % save at predefined time points
        if tInd == saveTimes(saveInd)
            pTag(saveInd, :)      = out(1:N);
            pTagSpine(saveInd, :) = out(N+1:end);
            saveInd               = saveInd + 1;
        end
    end
    clear done in out tInd

elseif strcmp(exp_imp, 'implicit')

    % coefficient matrix for implicit method
    % ---------------------------------------------------------------------
    
    % matrix size: (N+1 + N) x (N+1 + N)
    % add a fictitious node N+1 to represent the no flux boundary condition
    %
    %                   <-------- N+1 -------x------- N ------->
    %
    %  |   pSom   |     | 1 0  0             | 0               |     |  p(1, i+1)  |    ^
    %  | p(2, i)  |     | b a1 b             |   c             |     |  p(1, i+1)  |    |
    %  | p(3, i)  |     |   b a2 b           |     c           |     |  p(1, i+1)  |    |
    %  |   ...    |     |     .  . .         |       .         |     |     ...     |    |
    %  |   ...    |     |       .  . .       |         .       |     |     ...     |   N+1
    %  |   ...    |     |         .  . .     |           .     |     |     ...     |    |
    %  |   ...    |     |           .  . .   |             c   |     |     ...     |    |
    %  | p(N, i)  |     |             b aN b |               c |     |  p(N, i+1)  |    |
    %  |    0     |     |            -1 0  1 |               0 |     | p(N+1, i+1) |    v
    %  |          |  =  |--------------------------------------|  x  |             |  
    %  | pS(1, i) |     | d1                 | c               |     |  pS(1, i+1) |    ^
    %  | pS(1, i) |     |   d2               |   c             |     |  pS(1, i+1) |    |
    %  | pS(1, i) |     |     d3             |     c           |     |  pS(1, i+1) |    |
    %  |   ...    |     |        .           |       .         |     |     ...     |    |
    %  |   ...    |     |          .         |         .       |     |     ...     |    N
    %  |   ...    |     |            .       |           .     |     |     ...     |    |
    %  |   ...    |     |              .     |             c   |     |     ...     |    |
    %  | pS(N, i) |     |               dN 0 |               c |     |  pS(1, i+1) |    v
    
    %           2 D_p dT
    % aK = 1 + ---------- + lambda dT + supplyRatio(xK) uptakeRate dT
    %             dX^2
    %
    %         D dT
    % b  = - ------
    %         dX^2
    %  
    % c  = - exitRate dT
    %
    % d  = - supplyRatio(xK) uptakeRate dT
    %
    % e  = 1 + exitRate dT + lambda dT 
    
    % set up coefficient matrix for implicit method
    % ---------------------------------------------------------------------

    % left upper
    a            = 1 + 2*params.DEns_p*DeltaT/DeltaX^2 + params.lambda_p*DeltaT + supplyRatio.*params.uptakeRate.*DeltaT;
    b            = -params.DEns_p*DeltaT/DeltaX^2;
    b            = b.*ones(1, N);
    LU           = diag([b(2:N), 0], -1) + diag([1, a(2:N), 1], 0) + diag([0, b(2:N)], 1);
    LU(N+1, N-1) = -1;
    % right upper
    c            = -params.exitRate*DeltaT;
    c            = [0, c.*ones(1, N-1)];
    RU           = [diag(c, 0); zeros(1, N)];
    % left lower
    d            = -supplyRatio.*params.uptakeRate.*DeltaT;
    LL           = [diag(d, 0), zeros(N, 1)];
    % right lower
    e            = 1 + params.exitRate*DeltaT + params.lambda_p*DeltaT;
    e            = e*ones(1, N);
    RL           = diag(e, 0);
    % consolidate in one matrix
    A_imp        = sparse([LU, RU; LL, RL]);
    clear a b c d e LL LU RL RU

    % run implicit method
    % ---------------------------------------------------------------------

    % counter for progress of script
    done                          = 0;
    fprintf('\n %d percent done', done);

    in    = sparse(zeros(2*N+1, 1));
    in(1) = pSom;

    for tInd = 2:numel(t)
        % print progress 
        if round(tInd/numel(t)*100) > done
            done                  = round(tInd/numel(t)*100);
            fprintf('\n %d percent done', done);
        end
        % solve equation system
        out = A_imp\in;
        % re-set input array
        in      = out;
        in(N+1) = 0;
        % save at predefined time points
        if tInd == saveTimes(saveInd)
            pTag(saveInd, :)      = out(1:N);
            pTagSpine(saveInd, :) = out(N+2:end);
            saveInd               = saveInd + 1;
        end
    end
    clear done in out tInd
end

%% save results with a chosen date string

save(['files\', dateStr, '_Shank3Dynamics_results_', exp_imp], ...
     'params', 'x', 'pSom', 'pInit', 'pInitSpine', 'pTag', 'pTagSpine', 't', 'saveTimes')
