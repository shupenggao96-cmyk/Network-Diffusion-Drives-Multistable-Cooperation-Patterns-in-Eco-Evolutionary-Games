function dm_low = Func_Calculate_Critical_Dm(P)
%FUNC_CALCULATE_CRITICAL_DM Estimate the critical public-goods diffusivity.
%   DM_LOW = FUNC_CALCULATE_CRITICAL_DM(P) evaluates the lower threshold
%   for the public-goods diffusion coefficient dm based on the local
%   Jacobian and the continuous Laplacian-eigenvalue instability condition.
%
%   Input:
%       P      - Structure of model parameters, including P.du.
%
%   Output:
%       dm_low - Calculated diffusion threshold for dm.
%
%   This continuum-based threshold does not by itself establish a
%   Turing instability on a finite network: at least one nonzero graph
%   Laplacian eigenvalue must also fall within an unstable interval.

%% 1. Evaluate the equilibrium and local Jacobian
[ue, me] = Fucn_Get_Equilibrium(P);
[Fu, Fm, Gu, Gm] = Fucn_Get_Jacobian(P, ue, me);
du = P.du;

%% 2. Construct the discriminant condition for diffusion-driven instability
% The boundary between stable and unstable spectral intervals is given by:
%   (du*Gm + dm*Fu)^2 - 4*du*dm*(Fu*Gm - Fm*Gu) = 0.
%
% Rearranging yields A*dm^2 + B*dm + C = 0.
A = Fu^2;
B = -2 * du * Fu * Gm + 4 * du * Fm * Gu;
C = du^2 * Gm^2;

%% 3. Calculate candidate threshold values
% Find the two algebraic roots of the boundary condition.
p = [A, B, C];
roots_dm = roots(p);
dm_low_1 = max(roots_dm);

% The condition du*Gm + dm*Fu > 0 is also needed for an unstable
% interval at positive absolute Laplacian eigenvalues (when Fu > 0).
dm_low_2 = -du * Gm / Fu;

% Use the larger bound as the estimated onset threshold.
dm_low = max(dm_low_1, dm_low_2);

%% 4. Report the calculated threshold
fprintf('Critical diffusion threshold dm_low = %.4f (current du = %.4f)\n', ...
    dm_low, du);

end
