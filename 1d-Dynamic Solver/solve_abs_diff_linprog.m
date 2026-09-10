function [x_final, fval] = solve_abs_diff_linprog(D, A, b, lb, ub)
% Solves: min ||Dx||_1 subject to Ax = b, lb <= x <= ub
% Inputs:
%   D        - Difference matrix (m x n) where m is the number of differences
%   A, b     - Equality constraints (for x)
%   lb, ub   - Lower and upper bounds (for x)

[m, n] = size(D); % m is number of diff terms, n is length of x

% ---------------------------------------------------------
% 1. Set up the combined objective vector for [x; t]
% We want to minimize sum(t), so the weights on x are 0, and on t are 1.
% ---------------------------------------------------------
f_lin = [zeros(n, 1); ones(m, 1)];

% ---------------------------------------------------------
% 2. Set up the inequality constraints: Dx - t <= 0 and -Dx - t <= 0
% ---------------------------------------------------------
% Matrix multiplying [x; t]
A_ineq = [ D, -eye(m);
    -D, -eye(m)];

b_ineq = zeros(2*m, 1);

% ---------------------------------------------------------
% 3. Set up the equality constraints: Ax = b
% (t is not restricted by equality, so we pad A with zeros)
% ---------------------------------------------------------
A_eq = [A, zeros(size(A, 1), m)];
b_eq = b;

% ---------------------------------------------------------
% 4. Set up bounds for [x; t]
% t must be >= 0, and has no upper bound (Inf).
% ---------------------------------------------------------
lb_combined = [lb; zeros(m, 1)];
ub_combined = [ub; inf(m, 1)];

% ---------------------------------------------------------
% 5. Solve using linprog
% ---------------------------------------------------------

% options = optimoptions('linprog', ...
%     'Algorithm',        'dual-simplex-highs', ...
%     'Display',          'final', ...
%     'MaxIterations',    1e5, ...
%     'OptimalityTolerance', 1e-12, ...
%     'ConstraintTolerance', 1e-12);
% 
% fprintf('\nCalling linprog (dual-simplex-highs)...\n');
% [X_opt, fval, exitflag, output] = linprog(f_lin, A_ineq, b_ineq, A_eq, b_eq, lb_combined, ub_combined, options);
exitflag=0;
if exitflag < 1
    options = optimoptions('linprog', ...
        'Algorithm',        'interior-point', ...
        'Display',          'final', ...
        'MaxIterations',    1e4, ...
        'OptimalityTolerance', 1e-12, ...
        'ConstraintTolerance', 1e-12);

    fprintf('\nCalling linprog (interior-point)...\n');
    [X_opt, fval, exitflag, output] = linprog(f_lin, A_ineq, b_ineq, A_eq, b_eq, lb_combined, ub_combined, options);
end


% Extract the actual x variables from the combined solution
if exitflag == 1
    x_final = X_opt(1:n);
    % disp('Optimal solution found!');
else
    warning('linprog did not converge to an optimal solution. Check constraints.');
    x_final = [];
end
end