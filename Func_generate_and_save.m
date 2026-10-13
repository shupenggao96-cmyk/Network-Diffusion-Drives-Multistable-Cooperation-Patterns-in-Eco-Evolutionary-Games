function [deltau, deltam] = Func_generate_and_save(N_nodes, filename)
%FUNC_GENERATE_AND_SAVE Generate and save random initial perturbations.
%   [DELTAU, DELTAM] = FUNC_GENERATE_AND_SAVE(N_NODES, FILENAME)
%   generates independent standard-normal random perturbations for
%   cooperator density and public-goods level at each network node.
%
%   Inputs:
%       N_nodes  - Number of network nodes.
%       filename - Destination MAT-file path.
%
%   Outputs:
%       deltau - N_nodes-by-1 Gaussian perturbation vector for u.
%       deltam - N_nodes-by-1 Gaussian perturbation vector for m.
%
%   The random vectors are saved so that subsequent simulations can
%   reuse the same initial perturbation realization. A fixed random
%   seed is not set in this function.

%% 1. Generate independent standard-normal perturbations
deltau = randn(N_nodes, 1);
deltam = randn(N_nodes, 1);

%% 2. Save the perturbations for reuse in future simulations
save(filename, 'deltau', 'deltam');

end
