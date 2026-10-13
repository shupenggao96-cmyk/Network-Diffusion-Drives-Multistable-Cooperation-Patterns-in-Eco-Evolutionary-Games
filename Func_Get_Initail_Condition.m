function [u0, m0] = Func_Get_Initail_Condition(P, N_nodes)
% FUNC_GET_INITAIL_CONDITION Generate initial conditions for numerical simulations.
%
% This function initializes the cooperator density and public goods
% variables by adding small Gaussian random perturbations to the
% homogeneous equilibrium of the local reaction system.
%
% Inputs:
%   P       - Structure containing the model parameters.
%   N_nodes - Total number of nodes in the network.
%
% Outputs:
%   u0      - Initial cooperator density at each node (N_nodes x 1).
%   m0      - Initial public goods level at each node (N_nodes x 1).

%% 1. Calculate the homogeneous equilibrium
[ue, me] = Fucn_Get_Equilibrium(P);

%% 2. Generate perturbed initial conditions
% Add independent Gaussian perturbations with a standard deviation
% of 0.0005 to the homogeneous equilibrium at each network node.
perturbation_scale = 0.0005;

u0 = ue + perturbation_scale * randn(N_nodes, 1);
m0 = me + perturbation_scale * randn(N_nodes, 1);

end