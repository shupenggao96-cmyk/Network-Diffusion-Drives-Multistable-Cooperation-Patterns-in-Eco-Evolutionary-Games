function [Fu, Fm, Gu, Gm] = Fucn_Get_Jacobian(P, ue, me)
%FUCN_GET_JACOBIAN Evaluate the local reaction Jacobian at equilibrium.
%   [FU, FM, GU, GM] = FUCN_GET_JACOBIAN(P, UE, ME) calculates the
%   partial derivatives of the two local reaction functions F(u,m)
%   and G(u,m), evaluated at the homogeneous equilibrium (UE, ME).
%
%   Inputs:
%       P      - Structure of payoff and environmental parameters.
%       ue, me - Homogeneous equilibrium values of the two variables.
%
%   Outputs:
%       Fu - Partial derivative dF/du.
%       Fm - Partial derivative dF/dm.
%       Gu - Partial derivative dG/du.
%       Gm - Partial derivative dG/dm.
%
%   The Jacobian matrix is J = [Fu, Fm; Gu, Gm].

%% 1. Extract model parameters
% Payoff parameters.
a0 = P.a0; a1 = P.a1; b0 = P.b0; b1 = P.b1;
c0 = P.c0; c1 = P.c1; d0 = P.d0; d1 = P.d1;

% Environmental-feedback parameters and local population density.
alpha = P.alpha; beta = P.beta; delta = P.delta;
epsilon = P.epsilon; k = P.k; N = P.NA;

%% 2. Calculate coefficients used in the reaction functions
eta1 = a1 - a0 + b0 - b1 + c0 - c1 + d1 - d0;
eta2 = d0 - d1 + b1 - b0;
eta3 = b0 + c0 - a0 - d0;

%% 3. Evaluate the four Jacobian entries
% Derivatives of the cooperator-density reaction function F(u,m).
Fu = (b0 * delta * N * (N - 2 * ue) ...
    - d0 * delta * N * (N - 2 * ue) ...
    + eta1 * me * (2 * N - 3 * ue) * ue ...
    + eta2 * me * N * (N - 2 * ue) ...
    + delta * eta3 * ue * (-2 * N + 3 * ue)) / N^2;

Fm = ((N - ue) * ue * (eta2 * N + eta1 * ue)) / N^2;

% Derivatives of the public-goods reaction function G(u,m).
Gu = (alpha * beta * epsilon) / (alpha + ue)^2;
Gm = -epsilon * k * N;

end
