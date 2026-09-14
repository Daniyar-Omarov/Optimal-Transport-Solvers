function [cost, f_path, m_path, Pi_final] = PS_algoritm(Nx, Nt, mu1, mu2, xa, xb, ya, yb)
% Solves Benamou-Brenier Optimal Transport via Proximal Splitting (Asymmetric DR)
% Outputs:
%   cost     - Full W2^2 approximation
%   f_path   - [Nx, Nt+1] density along geodesic
%   m_path   - [Nx+1, Nt] flux along geodesic
%   Pi_final - [Nx, 1] terminal marginal
%   xhat     - [Nx, 1] Source midpoint grid
%   T_x      - [Nx, 1] Optimal map T(x) = x + v(x, t=0), evaluated at xhat

% 1. Determine Global Domain and Spatial Nodes for PDE Solver
x_min = min(xa, ya);
x_max = max(xb, yb);
L = x_max - x_min;

dx = L / Nx;
dt = 1.0 / Nt;
x_f = x_min + ((0.5 : Nx - 0.5) * dx)'; % Global spatial cell centers

% Evaluate and mask boundary distributions to their exact supports
mu1_vec = zeros(Nx, 1);
mask1 = (x_f >= xa) & (x_f <= xb);
mu1_vec(mask1) = mu1(x_f(mask1));

mu2_vec = zeros(Nx, 1);
mask2 = (x_f >= ya) & (x_f <= yb);
mu2_vec(mask2) = mu2(x_f(mask2));

% Normalize (mass must = 1)
if sum(mu1_vec) > 0, mu1_vec = mu1_vec / (sum(mu1_vec)); end
if sum(mu2_vec) > 0, mu2_vec = mu2_vec / (sum(mu2_vec)); end

% Dimensions and index mapping
Nf = Nx * (Nt + 1);
Nm = (Nx + 1) * Nt;
Nvars = Nf + Nm;
idx_f = @(i, j) i + (j - 1) * Nx;
idx_m = @(i, j) Nf + i + (j - 1) * (Nx + 1);

% ---------------------------------------------------------
% 2. Build Constraint Matrix C (for Proj_C)
% ---------------------------------------------------------
I = []; J = []; V = [];
b = zeros(Nx*Nt + 2*Nx + 2*Nt, 1);
eq_idx = 1;

% Boundary conditions
for i = 1:Nx
    I(end+1) = eq_idx; J(end+1) = idx_f(i, 1); V(end+1) = 1;
    b(eq_idx) = mu1_vec(i); eq_idx = eq_idx + 1;
    I(end+1) = eq_idx; J(end+1) = idx_f(i, Nt+1); V(end+1) = 1;
    b(eq_idx) = mu2_vec(i); eq_idx = eq_idx + 1;
end
for j = 1:Nt
    I(end+1) = eq_idx; J(end+1) = idx_m(1, j); V(end+1) = 1;
    b(eq_idx) = 0; eq_idx = eq_idx + 1;
    I(end+1) = eq_idx; J(end+1) = idx_m(Nx+1, j); V(end+1) = 1;
    b(eq_idx) = 0; eq_idx = eq_idx + 1;
end

% Continuity Equation (Divergence)
for j = 1:Nt
    for i = 1:Nx
        I(end+1) = eq_idx; J(end+1) = idx_f(i, j+1); V(end+1) = 1/dt;
        I(end+1) = eq_idx; J(end+1) = idx_f(i, j);   V(end+1) = -1/dt;
        I(end+1) = eq_idx; J(end+1) = idx_m(i+1, j); V(end+1) = 1/dx;
        I(end+1) = eq_idx; J(end+1) = idx_m(i, j);   V(end+1) = -1/dx;
        b(eq_idx) = 0; eq_idx = eq_idx + 1;
    end
end

C = sparse(I, J, V, eq_idx-1, Nvars);
C(end, :) = []; b(end) = [];
L_CCT = decomposition(C * C', 'chol');

% ---------------------------------------------------------
% 3. Build Midpoint Interpolation Operator K (for Proj_{C_{s,c}})
% ---------------------------------------------------------
I_K = []; J_K = []; V_K = []; row_idx = 1;
for j = 1:Nt
    for i = 1:Nx
        I_K(end+1) = row_idx; J_K(end+1) = idx_f(i, j); V_K(end+1) = 0.5;
        I_K(end+1) = row_idx; J_K(end+1) = idx_f(i, j+1); V_K(end+1) = 0.5;
        row_idx = row_idx + 1;
    end
end
for j = 1:Nt
    for i = 1:Nx
        I_K(end+1) = row_idx; J_K(end+1) = idx_m(i, j); V_K(end+1) = 0.5;
        I_K(end+1) = row_idx; J_K(end+1) = idx_m(i+1, j); V_K(end+1) = 0.5;
        row_idx = row_idx + 1;
    end
end
K = sparse(I_K, J_K, V_K, 2*Nx*Nt, Nvars);
L_IKK = decomposition(speye(Nvars) + K' * K, 'chol');

% ---------------------------------------------------------
% 4. Proximal Splitting (Asymmetric Douglas-Rachford)
% ---------------------------------------------------------
U = zeros(Nvars, 1);
V = zeros(2*Nx*Nt, 1);
w_U = U; w_V = V;

% gamma = 1.0; alpha = 1.0; max_iters = 1e3;
gamma = 1/75; alpha = 1.998; max_iters = 1e4;

for iter = 1:max_iters
    U_bar = 2 * U - w_U;
    V_bar = 2 * V - w_V;

    lambda = L_CCT \ (C * U_bar - b);
    U_prox = U_bar - C' * lambda;

    fc_bar = V_bar(1:Nx*Nt);
    mc_bar = V_bar(Nx*Nt+1:end);

    f_star = max(fc_bar, 0) + gamma + abs(mc_bar);
    for k = 1:10
        val = (f_star - fc_bar) .* (f_star + gamma).^2 - 0.5 * gamma * mc_bar.^2;
        der = (f_star + gamma).^2 + 2 * (f_star - fc_bar) .* (f_star + gamma);
        f_star = f_star - val ./ der;
    end

    fc_prox = zeros(size(fc_bar));
    mc_prox = zeros(size(mc_bar));
    mask = f_star > 0;
    fc_prox(mask) = f_star(mask);
    mc_prox(mask) = f_star(mask) .* mc_bar(mask) ./ (f_star(mask) + gamma);
    V_prox = [fc_prox; mc_prox];

    w_U = w_U + alpha * (U_prox - U);
    w_V = w_V + alpha * (V_prox - V);

    U = L_IKK \ (w_U + K' * w_V);
    V = K * U;
end

% ---------------------------------------------------------
% 5. Format Outputs & Compute Optimal Map T(x) on xhat
% ---------------------------------------------------------
f_path = reshape(U(1:Nf), [Nx, Nt+1]);
m_path = reshape(U(Nf+1:end), [Nx+1, Nt]);
Pi_final = f_path(:, end);

fc = V(1:Nx*Nt);
mc = V(Nx*Nt+1:end);
pos_f = fc > 1e-9;

cost_vec = zeros(size(fc));
cost_vec(pos_f) = mc(pos_f).^2 ./ (fc(pos_f));
cost = sum(cost_vec) * dt;

end