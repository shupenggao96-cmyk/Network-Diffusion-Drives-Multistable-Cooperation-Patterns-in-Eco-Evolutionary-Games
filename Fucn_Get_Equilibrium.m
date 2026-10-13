function [ue, me] = Fucn_Get_Equilibrium(P)
%FUCN_GET_EQUILIBRIUM Compute the nontrivial homogeneous equilibrium.
%   [UE, ME] = FUCN_GET_EQUILIBRIUM(P) calculates the equilibrium
%   cooperator density UE and public-goods level ME of the local
%   eco-evolutionary game model in the absence of diffusion.
%
%   Input:
%       P  - Structure containing payoff and environmental parameters.
%            P.NA is the fixed population density at each node.
%
%   Outputs:
%       ue - Equilibrium cooperator density.
%       me - Equilibrium public-goods level.
%
%   The function selects the branch with the minus sign in the
%   quadratic formula. The admissibility of the resulting equilibrium
%   (e.g., 0 < ue < N) is not checked here.

%% 1. Extract model parameters
% Payoff parameters.
a0 = P.a0; a1 = P.a1; b0 = P.b0; b1 = P.b1;
c0 = P.c0; c1 = P.c1; d0 = P.d0; d1 = P.d1;

% Environmental-feedback parameters and local population density.
alpha = P.alpha; beta = P.beta; delta = P.delta;
k = P.k; N = P.NA;

%% 2. Calculate coefficients of the local reaction dynamics
eta1 = a1 - a0 + b0 - b1 + c0 - c1 + d1 - d0;
eta2 = d0 - d1 + b1 - b0;
eta3 = b0 + c0 - a0 - d0;

%% 3. Construct the quadratic equilibrium equation
% The nontrivial equilibrium satisfies:
%   gamma1*ue^2 + gamma2*N*ue + gamma3*N^2 = 0.
gamma1 = beta * eta1 - N * delta * k * eta3;
gamma2 = ((b0 - d0) * N - eta3 * alpha) * k * delta + beta * eta2;
gamma3 = alpha * delta * k * (b0 - d0);

%% 4. Solve for the equilibrium cooperator density
% Evaluate the discriminant to check whether this quadratic has real roots.
discriminant = (gamma2 * N)^2 - 4 * gamma1 * (N^2 * gamma3);
if discriminant < 0
    error('The equilibrium equation has no real roots. Check the model parameters.');
end

% Select the root corresponding to the minus sign in the quadratic formula.
ue = (-gamma2 * N - sqrt(discriminant)) / (2 * gamma1);

%% 5. Recover the equilibrium public-goods level
% At equilibrium, G(ue, me) = 0 implies:
%   me = beta*ue / (k*N*(alpha + ue)).
me = (beta * ue) / (k * N * (alpha + ue));

end
