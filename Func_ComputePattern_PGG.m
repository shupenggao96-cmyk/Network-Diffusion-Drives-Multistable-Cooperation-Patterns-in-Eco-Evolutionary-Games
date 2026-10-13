function [U_stable, M_stable, U, M] = Func_ComputePattern_PGG(params, N_pop, L, T, deltaT, U0, M0)
%FUNC_COMPUTEPATTERN_PGG Simulate a networked eco-evolutionary public-goods game.
%   [U_STABLE, M_STABLE, U, M] = FUNC_COMPUTEPATTERN_PGG(...)
%   integrates the coupled reaction-diffusion equations using the
%   explicit forward Euler method.
%
%   Inputs:
%       params - Structure containing the payoff parameters,
%                environmental parameters, population density (NA),
%                and diffusion coefficients (du and dm).
%       N_pop  - Number of network nodes.
%       L      - Network Laplacian matrix of size N_pop-by-N_pop.
%                The diffusion terms use the convention L = A - D.
%       T      - Total integration time.
%       deltaT - Integration time step.
%       U0     - Initial cooperator densities (N_pop-by-1).
%       M0     - Initial public-goods levels (N_pop-by-1).
%
%   Outputs:
%       U_stable - Final cooperator densities at the end of integration.
%       M_stable - Final public-goods levels at the end of integration.
%       U        - Cooperator-density history; rows are nodes and
%                  columns are recorded time points.
%       M        - Public-goods history with the same layout as U.
%
%   The state is updated at every deltaT, but trajectories are stored
%   only at integer-valued simulation times to limit memory usage.
%   The final states are not necessarily equilibrium states.

%% 1. Extract model parameters
% Payoff parameters.
a0 = params.a0; a1 = params.a1; b0 = params.b0; b1 = params.b1;
c0 = params.c0; c1 = params.c1; d0 = params.d0; d1 = params.d1;

% Environmental feedback and network diffusion parameters.
delta = params.delta; eps = params.epsilon;
beta = params.beta; alpha = params.alpha; k = params.k;
du = params.du; dm = params.dm;
N = params.NA;  % Fixed population density at each network node.

% Coefficients appearing in the local cooperator-density dynamics.
eta1 = a1 - a0 + b0 - b1 + c0 - c1 + d1 - d0;
eta2 = d0 - d1 + b1 - b0;
eta3 = b0 + c0 - a0 - d0;

%% 2. Allocate storage and set initial conditions
% Construct the integration grid and the integer-time recording grid.
timeGrid = 0:deltaT:T;
record_steps = 0:1:T;

% Store each network node in a row and each recorded time in a column.
U = zeros(N_pop, length(record_steps));
M = zeros(N_pop, length(record_steps));

% Initialize the state vectors and store the initial state at t = 0.
u_curr = U0;
m_curr = M0;
U(:, 1) = u_curr;
M(:, 1) = m_curr;
record_idx = 1;

%% 3. Integrate the reaction-diffusion equations by forward Euler
for n = 2:length(timeGrid)
    % Evaluate the local reaction terms independently at each node.
    F_val = reaction_F(u_curr, m_curr, eta1, eta2, eta3, delta, d0, b0, N);
    G_val = reaction_G(u_curr, m_curr, eps, beta, alpha, k, N);

    % Apply network diffusion and advance both state variables.
    % With L = A - D, L*u_curr represents the diffusive coupling
    % between a node and its network neighbors.
    u_next = u_curr + (F_val + du * L * u_curr) * deltaT;
    m_next = m_curr + (G_val + dm * L * m_curr) * deltaT;

    % Update the state for the next integration step.
    u_curr = u_next;
    m_curr = m_next;

    % Record the state whenever the simulation time is an integer.
    if mod(timeGrid(n), 1) == 0
        record_idx = record_idx + 1;
        U(:, record_idx) = u_curr;
        M(:, record_idx) = m_curr;
    end
end

%% 4. Return the final state of the integration
U_stable = u_curr;
M_stable = m_curr;

%% 5. Define the local reaction functions
    function val = reaction_F(u, m, e1, e2, e3, delta, d0, b0, N)
        % Cooperator-density reaction term:
        % F(u,m) = (u/N^2)*(N-u)*[(e1*u + N*e2)*m
        %          - delta*(e3*u + N*(d0-b0))].
        term1 = (u ./ N^2) .* (N - u);
        term2 = (e1 .* u + N * e2) .* m ...
            - delta * (e3 .* u + N * (d0 - b0));
        val = term1 .* term2;
    end

    function val = reaction_G(u, m, eps, beta, alpha, k, N)
        % Public-goods reaction term:
        % G(u,m) = eps*[beta*u/(alpha+u) - k*N*m].
        val = eps * ((beta .* u) ./ (alpha + u) - k * N * m);
    end

end
