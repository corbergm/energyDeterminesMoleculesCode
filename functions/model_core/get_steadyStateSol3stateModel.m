function solution = get_steadyStateSol3stateModel(params) 

% author: CORNELIUS BERGMANN
%
% This function calculates the numerical solution of the following ODE
% system on the interval [0, L]:
%
% (EQ I)   0 = - vPlus  yP' - bPlus  yP + alpha y0
% (EQ II)  0 =   vMinus yM' - bMinus yM + alpha y0
% (EQ III) 0 =    bPlus yP  + bMinus yM + D y0'' - 2 alpha y0 - lambda y0
%
% which can be written as (using a := alpha)
%
% | yP' |   | -bPlus/vPlus,             0,   0,    a/vPlus |   | yP |
% | yM' |   |            0, bMinus/vMinus,   0,  -a/vMinus |   | yM |
% | y1' | = |       -bPlus,       -bMinus,   0, 2*a+lambda | x | y1 |
% | y0' |   |            0,             0, 1/D,          0 |   | y0 |
% 
% using y1 := D y0'. For the coefficient matrix we calculate the
% determinant to 
%
%        lambda     bPlus bMinus
% det = -------- x -------------- > 0.
%          D        vPlus vMinus 
% 
% At the boundaries we demand with "atSoma > 0"
%
% (BC I)   0 =       yP (0) +        yM (0) + y0(0) - atSoma
% (BC II)  0 = vPlus yP'(0) - vMinus yM'(0)
% (BC III) 0 = vPlus yP'(L) - vMinus yM'(L)
% (BC IV)  0 =       y1 (L) .
%
% For (BC II), (BC III) we can now write using (EQ I), (EQ II)
%
% (BC II)  0 = bPlus yP(0) - bMinus yM(0)
% (BC III) 0 = bPlus yP(L) - bMinus yM(L)
%
% Input: 
% params   - parameter structure containing the fields:
%            "vPlus", "vMinus", "D", "length", "bPlus", "bMinus", "alpha",
%            "transpProb", "halflife", "atSoma", "maxRuns"
%
% Output: 
% solution - structure containing the solution of the ODE system with
%            additional fields "nMeshPts" (number of mesh points) and 
%            "maxErr" (log10 solver error)
%
% Note: The dendrite length "length" has to be smaller than or a multiple
%       of 100, i.e., 1-100, 200, 300, 400, ... are possible values.

%% Set model parameters

vPlus  = params.vPlus;
vMinus = params.vMinus;
D      = params.D;
L      = params.length;
bPlus  = params.bPlus;
bMinus = params.bMinus;
a      = params.alpha;
trProb = params.transpProb;
lambda = log(2)/params.halflife;
atSoma = params.atSoma;

%% Solve the ODE's

% Solve the ODE system given by 'transporteq' with boundary conditions 
% 'transportbc' on [0, L] directly. 
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
% Create and initiate the guess on a 100 point mesh on [0, L] 
xInit             = linspace(0, L, 100);
yInit             = guess(xInit);
solinit           = bvpinit(xInit, yInit);
% solve equations on [0, L] 
solution          = get_odeSolution(solinit, L, 0, 0);
% store final number of mesh points and solver error
solution.nMeshPts = solution.stats.nmeshpoints;
solution.maxErr   = log10(solution.stats.maxerr);

%% Local functions 

% ----------------------------- Initial guess -----------------------------

    function y = guess(x)
        % We use the following constant initial guess:
        %
        %  yP  = transpProb/2 x atSoma
        %  yM  = transpProb/2 x atSoma
        %  y0' = 0
        %  y0  = (1-transpProb) x atSoma
        %
        y = [(trProb/2)*atSoma
             (trProb/2)*atSoma
             0
             (1 - trProb)*atSoma];
    end

% ------------------------------- Equations -------------------------------
    
    function dydt = transporteq(t,y)
        dydt = [-bPlus/vPlus,             0,   0,    a/vPlus; ...
                           0, bMinus/vMinus,   0,  -a/vMinus; ...
                      -bPlus,       -bMinus,   0, 2*a+lambda; ...
                           0,             0, 1/D,          0] ...
                 *[y(1); y(2); y(3); y(4)];
    end

% -------------------------- Boundary conditions --------------------------

    function res = transportbc(ya, yb)
        res = [ya(1) + ya(2) + ya(4) - atSoma;
               ya(1) - ya(2)    ;
               yb(1) - yb(2)    ;
               yb(3)                         ];
    end

% -------------- Recursive call of bvp5c to increase accuracy -------------

    function solution = get_odeSolution(solinit, length, counter, nIter)
        % This function solves the ODE system given by 'transporteq' with
        % boundary conditions 'transportbc'. It uses the solution as the new
        % initial guess until either 
        %  1  the required tolerance is met or
        %  2  the maximal number of iterations 'counter' has been used
        
        % counter of iterations performed with given tolerance
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