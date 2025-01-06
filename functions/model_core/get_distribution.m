function solution = get_distribution(params) 

% author: CORNELIUS BERGMANN
%
% modified 22.01.2023
%
% This function calculates the distribution of mRNA and protein on a linear
% (without branching points) dendritic segment of length L:
%
%  (EQ 1)  0 = D_m m'' - lambda_m m
%  (EQ 2)  0 = D_p p'' - lambda_p p + tau m - uptake(p)
%
% using a continuous uptake function for protein
%
%                               perm p
%  uptake(p) = lambda_p -----------------------
%                         /  perm   \
%                        | --------- | x p + 1
%                         \ rho eta /
%
% with permeability perm = u_p/(nu_p + lambda_p), where u_p is the 
% switching rate of protein from the dendrite into the spine and
% nu_p is the switching rate of the inverse process. 
% Boundary conditions are 
%
%                           1       phi
%  (BC 1)      0 = p(L) - ------ --------- eta_p ,
%                          perm   1 - phi
%             _
%            /                    somaRet     tau D_m
%           |  0 = D_p p'(0) - ------------- ---------- m'(0), somaRet < 1
%  (BC 2)  <                    1 - somaRet   lambda_m
%           |  0 = m'(0)                                     , somaRet = 1
%            \_
%
%  (BC 3)      0 = m'(L) ,
%  (BC 4)      0 = p'(L) ,
%
% where 0 < phi < 1 is the supply ratio of spines. 
% The calculation is straightforward and the solver is bvp5c.
%
% Input: 
% params   - parameter structure containing the fields:
%            "length", "DEns_m", "DEns_p", "tRate", "halflife_m",
%            "halflife_p", "permeability", "rhoSpines", "eta_p", "phi",
%            "somaRet", "maxRuns"
%
% Output: 
% solution - structure containing the solution of the ODE system with
%            additional fields "nMeshPts" (number of mesh points) and 
%            "maxErr" (log10 solver error)
%
% Note: The dendrite length "length" has to be smaller than or a multiple
%       of 250, i.e., 1-250, 500, 750, 1000, ... are possible values.

%% Set spine and particle parameters  

% load necessary parameters
L          = params.length;
D_m        = params.DEns_m;
D_p        = params.DEns_p;
tRate      = params.tRate;
lambda_m   = log(2)/params.halflife_m;
lambda_p   = log(2)/params.halflife_p;
perm       = params.permeability;
rho        = params.rhoSpines;
eta_p      = params.eta_p;
phi        = params.phi;
somaRet    = params.somaRet;

%% Solve the ODE's   

% Solve the ODE system given by 'transporteq' with boundary conditions 
% 'transportbc' for increasing dendrite lengths. The equations are 
% iteratively solved on 250, 500, 750, 1000 micron segments up to the 
% maximal length 'L'. Solutions obtained on one interval are then extended 
% to the next longer interval till 'L' is reached. If 'L' < 250, the ODE
% system is solved on [0, L] directly. 
% Initially, 100 mesh points are used. bvp5c allows for NMax/4 mesh points 
% within [0, L] for a system with four equations.

% bvp5c is a finite difference code that implements the four-stage Lobatto
% IIIa formula. This is a collocation formula and the collocation 
% polynomial provides a C1-continuous solution that is fifth-order accurate 
% uniformly in [a,b]. The formula is implemented as an implicit Runge-Kutta 
% formula.

% Set the tolerances of bvp5c such that ||res/solution - e-8|| < e-8, for
% further definitions of absolute and relative tolerance see mathworks.com
% Also set maximal total number of mesh points to 2000.
options           = bvpset('AbsTol', 10^(-8), ...
                           'RelTol', 10^(-4), ...
                           'Stats' , 'off', ...
                           'NMax'  , 2000);
% Create and initiate the guess on a 100 point mesh on [0, min(L, 250)] 
LInit             = min([L, 250]);
xInit             = linspace(0, LInit, 100);
yInit             = guess(xInit);
solinit           = bvpinit(xInit, yInit);
% solve equations initially on a 'LInit' microns long dendrite
tmpL              = LInit;
solution          = get_odeSolution(solinit, tmpL, 0, 0);
% increase length and extend solution to longer dendrite
while tmpL < L    
    tmpL          = min([L, tmpL + 250]);
    solution      = get_odeSolution(bvpinit(solution, [0, tmpL]), tmpL, 0, 0);
end
% store final number of mesh points and solver error
solution.nMeshPts = solution.stats.nmeshpoints;
solution.maxErr   = log10(solution.stats.maxerr);

%% Local functions 

% ----------------------------- Initial guess -----------------------------

    function y = guess(x)
        % We use the following constant initial guess:
        %
        %        / 1,   if   somaRet = 0,
        %  m' = <
        %        \ 0,   if   somaRet = 1,
        %
        %         /   p         1      1       phi
        %        | ------- = ------- ------ --------- eta_p,   if   somaRet = 0,
        %  m  = <   tRate     tRate   perm   1 - phi
        %        |  
        %         \                                       0,   if   somaRet = 1,
        %
        %  p' = 1,
        %
        %         1       phi
        %  p  = ------ --------- eta_p ,
        %        perm   1 - phi
        %  
        y = [(1 - somaRet)
             (1 - somaRet).*phi./(perm.*(1 - phi)).*eta_p./tRate
             1
             phi./(perm.*(1 - phi)).*eta_p];
    end

% ------------------------------- Equations -------------------------------

    function dydt = transporteq(t, y)
        % Starting with the ODE system 
        %
        %   0 = D_m m'' - lambda_m m 
        %   0 = D_p p'' - lambda_p p + tRate m - uptake(p)
        %
        % on [0, L]. Define
        %
        %               f = D_m m'
        %               g = D_p p'
        %
        % what leads to the system of first order ODE's
        %
        %   f' =         lambda_m m
        %   m' = 1/D_m f
        %   g' =           -tRate m         lambda_p p  + uptake(p)
        %   p' =                    1/D_p g
        %
        % with
        %                               perm p
        %  uptake(p) = lambda_p -----------------------  .
        %                         /  perm   \
        %                        | --------- | x p + 1
        %                         \ rho eta /
        %
        % In matrix notation:
        dydt = [    0, lambda_m,     0,        0; ...
                1/D_m,        0,     0,        0; ...
                    0,   -tRate,     0, lambda_p; ...
                    0,        0, 1/D_p,        0] ...
              *[y(1); y(2); y(3); y(4)];
        % Add protein uptake:
        dydt = dydt + [0; 0; lambda_p*perm*y(4)/(perm/(rho*eta_p)*y(4) + 1); 0];
    end

% -------------------------- Boundary conditions --------------------------

    function res = transportbc(ya, yb)
        %
        %                           1       phi
        %  (BC 1)      0 = p(L) - ------ --------- eta_p ,
        %                          perm   1 - phi
        %             _
        %            /               somaRet       tau 
        %           |  0 = g(0) - ------------- ---------- f(0), somaRet < 1
        %  (BC 2)  <               1 - somaRet   lambda_m
        %           |  0 = f(0)                                , somaRet = 1
        %            \_
        %
        %  (BC 3)      0 = f(L) ,
        %  (BC 4)      0 = g(L) ,
        %
        if somaRet < 1
            res = [yb(4) - phi/(perm*(1 - phi))*eta_p                ;
                   ya(3) - somaRet/(1 - somaRet)*tRate/lambda_m*ya(1);
                   yb(1)                                             ;
                   yb(3)                                             ];
        elseif somaRet == 1
            res = [yb(4) - phi/(perm*(1 - phi))*eta_p;
                   ya(1)                             ;
                   yb(1)                             ;
                   yb(3)                             ];
        end
    end

% -------------- Recursive call of bvp5c to increase accuracy -------------

    function solution = get_odeSolution(solinit, length, counter, nIter)
        % This function solves the ODE system given by 'transporteq' with
        % boundary conditions 'transportbc'. It uses the solution as the new
        % initial guess until either 
        %  1  the required tolerance is met or
        %  2  the maximal number of iterations 'counter' has been used
        
        % counter of iterations with given tolerance
        counter          = counter + 1;
        % overall number of performed iterations (with all tolerances used)
        nIter            = nIter + 1;
        % solve equations using the transferred initial guess
        solution         = bvp5c(@transporteq, @transportbc, solinit, options);
        % store "counter" and "nIter" within "solution" structure
        solution.counter = counter;
        solution.nIter   = nIter;
        % check, if actual solution meets the requested error tolerance 
        % which is set to 10 times the targeted "relative tolerance" of the 
        % solver bvp5c.
        % If the tolerance is not met after "maxRuns" iterations, reduce
        % the requested tolerance by a factor of 10 and reset the iteration 
        % counter.
        if solution.stats.maxerr > 10*options.RelTol
            if counter < params.maxRuns
                solution = get_odeSolution(bvpinit(solution, [0, length]), length, counter, nIter);
            else
                options.RelTol = options.RelTol.*10;
                solution = get_odeSolution(bvpinit(solution, [0, length]), length, 0, nIter);
            end
        end
    end

end

