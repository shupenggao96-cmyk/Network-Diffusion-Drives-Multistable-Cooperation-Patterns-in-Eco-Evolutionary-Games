function is_stable = Func_Check_Turing_Preconditions(P)
%FUNC_CHECK_TURING_PRECONDITIONS Check local stability without diffusion.
%   IS_STABLE = FUNC_CHECK_TURING_PRECONDITIONS(P) tests whether the
%   homogeneous equilibrium is linearly stable in the reaction-only
%   system, a prerequisite for classical diffusion-driven Turing instability.
%
%   Input:
%       P         - Structure containing the local model parameters.
%
%   Output:
%       is_stable - Logical true if trace(J) < 0 and det(J) > 0;
%                   logical false otherwise.
%
%   This function does not test the graph-specific diffusion-driven
%   instability condition or verify the equilibrium's physical range.

%% 1. Evaluate the equilibrium and its Jacobian
[ue, me] = Fucn_Get_Equilibrium(P);
[Fu, Fm, Gu, Gm] = Fucn_Get_Jacobian(P, ue, me);

%% 2. Calculate the trace and determinant
% For a two-dimensional reaction system, linear asymptotic stability
% requires trace(J) < 0 and det(J) > 0.
Condition1 = Fu + Gm;           % Trace of the local Jacobian.
Condition2 = Fu * Gm - Fm * Gu; % Determinant of the local Jacobian.

%% 3. Display the Jacobian and stability conditions
fprintf('--- Jacobian matrix ---\n');
fprintf('Fu = %.4f, Fm = %.4f\n', Fu, Fm);
fprintf('Gu = %.4f, Gm = %.4f\n', Gu, Gm);

fprintf('--- Stability without diffusion ---\n');
fprintf('u_star = %.4f, m_star = %.4f\n', ue, me);
fprintf('Tr(J)  = %.4f (required: < 0)\n', Condition1);
fprintf('Det(J) = %.4f (required: > 0)\n', Condition2);

%% 4. Return the local stability result
if Condition1 < 0 && Condition2 > 0
    is_stable = true;
    fprintf(['Result: The equilibrium is locally stable without diffusion; ', ...
        'a Turing instability may be possible.\n']);
else
    is_stable = false;
    fprintf(['Result: The equilibrium does not satisfy the local ', ...
        'stability prerequisite for classical Turing instability.\n']);
end

end
