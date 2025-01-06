function y = get_analyticalSteadyStateSol1StateModel(D, x, lambda, L, y0)

% author: CORNELIUS BERGMANN
%
% Input: 
% D        - diffusion constant
% x        - x-coordinates to evaluate 1-state model on
% lambda   - decay constant log(2)/T_{1/2}
% L        - interval length
% y0       - concentration at x0 = 0
%
% Output: 
% y        - analytical solution of the 1-state ODE system
%
% ---------------------------------------------------------------------
%
% This function implements the analytical solution of the following ODE
% on the interval [0, L]:
%
% (EQ I)   0 = D y'' - lambda y
%
% with the boundary conditions
%
% (BC I)   0 = y'(L)
% (BC II)  0 = y(0) - y0
%
% The solution is then given by:
%
%           __________             __________             __________
%  y0 exp(-V lambda/D  x) [ exp(2 V lambda/D  L) + exp(2 V lambda/D  x) ]
% ------------------------------------------------------------------------
%                          __________       
%                   exp(2 V lambda/D  L ) + 1

y = y0.*( exp(-sqrt(lambda/D)*x ).*( exp(2*sqrt(lambda/D)*L) + exp(2*sqrt(lambda/D).*x) ) ./ ( exp(2*sqrt(lambda/D)*L) + 1 ) );
