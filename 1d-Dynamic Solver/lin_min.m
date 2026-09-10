function [cost, Pi_final, Pi, diag] = lin_min(Nt, Nx, Ny, xa, xb, ya, yb)
% LIN_MIN  Solve the dynamic OT problem (flux formulation, quadratic cost).
%
% Inputs:
%   Nt  - number of time points (time grid: t_k=(k-1)*dt, k=1..Nt)
%   Nx  - number of z1 cells  (face grid: x_i=(i-1)*dx, i=1..Nx+1)
%   Ny  - number of z2 cells  (face grid: y_j=(j-1)*dy, j=1..Ny+1)
%
% Solution vector layout: sol = [Gamma; M1; M2]  (x->y->t, column-major)
%
%   Gamma_{i,j,k}: Nx*Ny*Nt   variables
%   M1_{i,j,k}:   (Nx+1)*Ny*(Nt-1) variables
%   M2_{i,j,k}:   Nx*(Ny+1)*(Nt-1) variables

%% Grid
x  = linspace(xa,xb,Nx+1)';          % z1 face points, (Nx+1)×1
y  = linspace(ya,yb,Ny+1)';          % z2 face points, (Ny+1)×1
dx = x(2)-x(1);  dy = y(2)-y(1);  dt = 1/(Nt-1);

xhat = 0.5*(x(1:end-1)+x(2:end));  % z1 cell midpoints, Nx×1
yhat = 0.5*(y(1:end-1)+y(2:end));  % z2 cell midpoints, Ny×1

% total_matches = sum(ismember(xhat, yhat)); [size(xhat,1) size(yhat,1) total_matches]

%% Marginal densities (evaluated at cell midpoints, then normalised)
mu1 = arrayfun(@mu1_func, xhat);  mu1 = mu1/sum(mu1);   % Nx×1
mu2 = arrayfun(@mu2_func, yhat);  mu2 = mu2/sum(mu2);   % Ny×1

%% Variable sizes
NGamma = Nx*Ny*Nt;
NM1    = (Nx+1)*Ny*(Nt-1);
NM2    = Nx*(Ny+1)*(Nt-1);
Ntot   = NGamma + NM1 + NM2;

%% Constraints and objective matrix
[Aeq, beq] = constraints(x, y, mu1, mu2, Nt, Nx, Ny);

D_mat = sparse(D_mat_gen(Nt, Nx, Ny));   % (Nx*Ny*(Nt-1)) × Ntot

%% Weight vector s_{i,j} = 2*|xhat_i - yhat_j|, reshaped to match D_mat rows
% D_mat row (i,j,k) <-> index i + (j-1)*Nx + (k-1)*Nx*Ny
s_ij = zeros(Nx, Ny);
for i = 1:Nx
    for j = 1:Ny
        s_ij(i,j) = 2*abs(xhat(i) - yhat(j));
    end
end
% Tile over k=1..Nt-1 (same spatial weights at every time step)
s_vec = repmat(s_ij(:), Nt-1, 1);    % (Nx*Ny*(Nt-1))×1

%% Bounds: Gamma >= 0; M1, M2 unconstrained
lb           = -Inf*ones(Ntot,1);
lb(1:NGamma) = 0;
ub           = Inf*ones(Ntot,1);

%% Solve via LP reformulation (min sum(t) s.t. s.*D*sol - t <= 0, etc.)
% Pass s_vec into solve_abs_diff_linprog so the LP correctly minimises
%   sum_k sum_{i,j} s_{i,j} * |(avg_M1 - avg_M2)_{i,j,k}|
% by using (s_vec.*D_mat) as the weighted difference matrix.
D_weighted = spdiags(s_vec,0,length(s_vec),length(s_vec)) * D_mat;

[sol, fval] = solve_abs_diff_linprog(D_weighted, Aeq, beq, lb, ub);

%% Evaluate objective and extract solution arrays
const_error = Aeq*sol - beq;

% Objective: sum of s_{i,j}*|avg_M1-avg_M2|*dx*dy*dt
cost    = zeros(2,1);
cost(1) = sum(abs(D_weighted*sol)) * dt;

% Reshape solution (x->y->t column-major: reshape to [x,y,t])
Pi   = reshape(sol(1:NGamma),                  Nx, Ny, Nt);
M1   = reshape(sol(NGamma+1       : NGamma+NM1),    Nx+1, Ny,   Nt-1);
M2   = reshape(sol(NGamma+NM1+1   : end),           Nx,   Ny+1, Nt-1);

Pi_final = Pi(:,:,end);   % Nx × Ny

% Alternative cost: W2^2 from the terminal coupling directly
dist_sq = (xhat - yhat').^2;           % Nx×Ny matrix of (xhat_i-yhat_j)^2
cost(2)  = sum(sum(Pi_final .* dist_sq));

%% Diagnostics
diag.const_err_max   = max(abs(const_error));
diag.row_sum_err_max = max(abs(sum(Pi_final,2) - mu1));   % z1-marginal
diag.col_sum_err_max = max(abs(sum(Pi_final,1)' - mu2));  % z2-marginal

end