%% Networked Eco-Evolutionary Game: Pattern Formation and Evolution
% This script investigates the spatiotemporal dynamics of cooperation
% in a networked eco-evolutionary public goods game.
%
% The simulation consists of two consecutive stages:
%   Stage 1: Evolution from a weakly perturbed homogeneous equilibrium.
%   Stage 2: Continued evolution after introducing additional
%            perturbations to the final state of Stage 1.
%
% Both stages use the same model parameters and diffusion coefficients.

clear; clc;
close all;

%% 1. Load the network Laplacian matrix
network_type = 'BA';       % Network topology (Barabasi-Albert network)
average_degree = 10;       % Average node degree
N_nodes = 100;             % Total number of network nodes

% Construct the path to the pre-generated network data.
network_root = ['./generateNetwork/LaplacianMatrixSort_', ...
    num2str(N_nodes), '/', ...
    network_type, num2str(average_degree), '.mat'];

% Load the network Laplacian matrix L.
load(network_root);

%% 2. Define model parameters
% Payoff parameters of the eco-evolutionary public goods game.
P.a0 = 2.9;  P.a1 = 3;
P.b0 = 2.3;  P.b1 = 1;
P.c0 = 0.1;  P.c1 = 2;
P.d0 = 1;    P.d1 = 2;

% Parameters governing environmental feedback and resource dynamics.
P.alpha = 0.9;
P.beta = 2.5;
P.delta = 2;
P.epsilon = 1.5;
P.k = 0.2;

% Fixed total population density at each network node.
P.NA = 3;

% Diffusion coefficients for cooperators and public goods.
P.du = 0.03;
P.dm = 0.85;

%% 3. Analyze stability and Turing instability conditions
% Check the local stability of the homogeneous equilibrium
% in the absence of diffusion.
is_stable = Func_Check_Turing_Preconditions(P);

% Calculate the critical public goods diffusion coefficient
% associated with diffusion-driven instability.
dm_critical = Func_Calculate_Critical_Dm(P);

% Compute the theoretical dispersion relation and the growth rates
% associated with the discrete network Laplacian eigenvalues.
[lambda_fine, sigma_max_fine, lambda_abs, sigma_net] = ...
    Func_Plot_Network_Turing_Dispersion(P, L);

%% 4. Stage 1: Evolution from a perturbed homogeneous equilibrium
% Initialize the system near the homogeneous equilibrium using
% small Gaussian random perturbations.
[u0, m0] = Func_Get_Initail_Condition(P, N_nodes);

% Set the simulation duration and numerical integration time step.
T = 399;
deltaT = 0.01;

% Solve the networked reaction-diffusion system using
% the explicit Forward Euler method.
[U_final, M_final, U_history, M_history] = ...
    Func_ComputePattern_PGG(P, N_nodes, L, T, deltaT, u0, m0);

%% 5. Stage 2: Evolution under additional perturbations
% Continue the simulation from the final state of Stage 1.
% Additional independent uniform random perturbations in [0, 0.1]
% are introduced to both state variables at each network node.
%
% The model parameters and diffusion coefficients remain unchanged.

T = 299;

% Construct the perturbed initial conditions for Stage 2.
u0 = U_final + 0.1 * rand(N_nodes, 1);
m0 = M_final + 0.1 * rand(N_nodes, 1);

% Continue the numerical integration with the perturbed states.
[U_final2, M_final2, U_history2, M_history2] = ...
    Func_ComputePattern_PGG(P, N_nodes, L, T, deltaT, u0, m0);

%% 6. Combine the simulation results
% Concatenate the trajectories from both simulation stages.
% Each row corresponds to a recorded time point, while each
% column represents an individual network node.
data1 = [U_history, U_history2]';
data2 = [M_history, M_history2]';

%% 7. Visualize the temporal evolution of cooperator density
% Plot the cooperator density trajectories across all network nodes.
% The horizontal axis represents the recorded simulation index.

figure;
plot(1:size(data1, 1), data1);

xlabel('Recorded time index');
ylabel('Cooperator density');
title('Temporal evolution of cooperator density');