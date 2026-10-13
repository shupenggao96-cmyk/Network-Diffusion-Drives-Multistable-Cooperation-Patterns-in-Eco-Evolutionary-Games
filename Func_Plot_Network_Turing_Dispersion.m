function [lambda_fine, sigma_max_fine, lambda_abs, sigma_net, lambda_cross] = ...
    Func_Plot_Network_Turing_Dispersion(P, L)
%FUNC_PLOT_NETWORK_TURING_DISPERSION Evaluate and plot the dispersion relation.
%   [LAMBDA_FINE, SIGMA_MAX_FINE, LAMBDA_ABS, SIGMA_NET, LAMBDA_CROSS] =
%   FUNC_PLOT_NETWORK_TURING_DISPERSION(P, L) calculates the maximum
%   linear growth rate for a continuous range of Laplacian eigenvalue
%   magnitudes and for the discrete eigenvalues of the supplied network.
%
%   Inputs:
%       P              - Model parameters, including diffusion coefficients
%                        P.du and P.dm.
%       L              - Network Laplacian matrix, using L = A - D.
%
%   Outputs:
%       lambda_fine    - Sampled nonnegative eigenvalue magnitudes.
%       sigma_max_fine - Maximum real growth rate on lambda_fine.
%       lambda_abs     - Sorted magnitudes of graph Laplacian eigenvalues.
%       sigma_net      - Maximum real growth rates at graph eigenvalues.
%       lambda_cross   - Detected intersections where sigma(lambda) = 0.
%
%   The continuous dispersion curve represents a spectral envelope.
%   Actual instability on a finite graph requires an unstable graph
%   eigenmode. Crossings are sought only within the plotted range.

%% 1. Evaluate the local Jacobian at the homogeneous equilibrium
du = P.du;
dm = P.dm;

[ue, me] = Fucn_Get_Equilibrium(P);
[Fu, Fm, Gu, Gm] = Fucn_Get_Jacobian(P, ue, me);

%% 2. Compute the network Laplacian eigenvalues
% Under L = A - D, the eigenvalues of an undirected graph are
% nonpositive. Their magnitudes are plotted on a nonnegative lambda axis.
Lambda_net = eig(full(L));
lambda_abs = sort(abs(real(Lambda_net)));

%% 3. Define the maximum linear growth-rate function
% For each nonnegative lambda, the linearized mode matrix is:
%   J(lambda) = [Fu - du*lambda, Fm;
%                Gu, Gm - dm*lambda].
% The larger real part of its two eigenvalues defines sigma(lambda).
sigma_function = @(l) real( ...
    ((Fu - du*l) + (Gm - dm*l) + ...
    sqrt(complex( ...
    ((Fu - du*l) + (Gm - dm*l)).^2 ...
    - 4 * ((Fu - du*l).*(Gm - dm*l) - Fm*Gu)))) / 2);

%% 4. Sample the continuous dispersion relation
% Extend the sampled interval beyond the largest network eigenvalue.
lambda_max = max(lambda_abs) * 1.2;
lambda_fine = linspace(0, lambda_max, 2000);
sigma_max_fine = arrayfun(sigma_function, lambda_fine);

%% 5. Locate zero crossings of the dispersion relation
lambda_cross = [];

% Identify intervals in which successive samples have opposite signs.
cross_index = find( ...
    sigma_max_fine(1:end-1) .* sigma_max_fine(2:end) < 0);

for i = 1:length(cross_index)
    index = cross_index(i);
    lambda_left = lambda_fine(index);
    lambda_right = lambda_fine(index + 1);

    % Refine each bracketed zero crossing with a scalar root solver.
    lambda_root = fzero(sigma_function, [lambda_left, lambda_right]);
    lambda_cross = [lambda_cross; lambda_root];
end

% Also retain sampled points that are numerically close to zero.
zero_tolerance = 1e-8;
zero_index = find(abs(sigma_max_fine) < zero_tolerance);

if ~isempty(zero_index)
    lambda_cross = [lambda_cross; lambda_fine(zero_index(:))];
end

% Merge numerical duplicates and sort the identified crossings.
if ~isempty(lambda_cross)
    lambda_cross = uniquetol(lambda_cross, 1e-6);
    lambda_cross = sort(lambda_cross);
end

%% 6. Evaluate growth rates at the discrete network eigenvalues
sigma_net = arrayfun(sigma_function, lambda_abs);

%% 7. Report the detected zero crossings
if isempty(lambda_cross)
    fprintf('No intersections with sigma = 0 were found in the sampled range.\n');
else
    fprintf('Intersections of sigma(lambda) with the lambda axis:\n');

    for i = 1:length(lambda_cross)
        fprintf('lambda_%d = %.8f\n', i, lambda_cross(i));
    end

    if length(lambda_cross) == 2
        fprintf('Continuous unstable interval: %.8f < lambda < %.8f\n', ...
            lambda_cross(1), lambda_cross(2));
    end
end

%% 8. Plot the continuous curve and discrete network modes
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;
grid on;
box on;

% Plot the spectral-envelope dispersion curve.
plot(lambda_fine, sigma_max_fine, ...
    'b-', 'LineWidth', 2, ...
    'DisplayName', 'Theoretical Curve');

% Overlay the growth rates corresponding to graph eigenvalues.
scatter(lambda_abs, sigma_net, 35, 'r', ...
    'filled', 'DisplayName', 'Network Eigenvalues');

% Draw the sigma = 0 reference line without using yline.
x_range = xlim;
plot(x_range, [0, 0], 'k--', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
xlim(x_range);

% Highlight and label the calculated intersections.
if ~isempty(lambda_cross)
    scatter(lambda_cross, zeros(size(lambda_cross)), ...
        80, 'ko', 'filled', ...
        'DisplayName', '\sigma(\lambda)=0');

    for i = 1:length(lambda_cross)
        text(lambda_cross(i), 0, ...
            sprintf('  \\lambda_{%d}=%.4f', i, lambda_cross(i)), ...
            'FontSize', 11, ...
            'Color', 'k', ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left');
    end
end

%% 9. Highlight growing network eigenmodes
% A positive growth rate indicates linear instability of that mode.
unstable_pts = sigma_net > 0;

if any(unstable_pts)
    plot(lambda_abs(unstable_pts), sigma_net(unstable_pts), ...
        'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
        'DisplayName', 'Unstable Modes');

    title('Turing Pattern Possible: Yes', 'Color', 'r');
else
    title(['Turing Pattern Possible: No ', ...
        '(No eigenvalues in unstable region)'], ...
        'Color', 'k');
end

xlabel('Laplacian Eigenvalue |\lambda_i|');
ylabel('Growth Factor Re(\sigma)');
legend('Location', 'best');

hold off;
end
