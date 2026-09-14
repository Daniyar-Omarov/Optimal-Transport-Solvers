function [cost, Pi_final, Pi, diag] = lin_min_2d(Nt, Nx, xa, xb)
% LIN_MIN_2D  Solve the dynamic OT problem (flux formulation, quadratic
%             cost) for two marginals mu1, mu2 on the square [xa,xb]^2,
%             i.e. z1, z2 each range over R^2.
%
% Inputs:
%   Nt        - number of time points (time grid: t_k=(k-1)*dt, k=1..Nt)
%   Nx        - number of cells PER SPATIAL DIRECTION, shared by all four
%               coordinates z1^1, z1^2, z2^1, z2^2 (both densities live on
%               the same square domain [xa,xb]^2)
%   xa, xb    - domain endpoints, shared by both densities (use 0,1 for
%               the unit square)
%   mu1_func  - function handle mu1_func(x,y) giving the (unnormalised)
%               density of mu1 at the point (x,y)
%   mu2_func  - function handle mu2_func(x,y) giving the (unnormalised)
%               density of mu2 at the point (x,y)
%
% Solution vector layout: sol = [Gamma; M11; M12; M21; M22]
%   Gamma_{i,j,a,b,k}: Nx^4*Nt         variables
%   M11_{i,j,a,b,k}:  (Nx+1)*Nx^3*(Nt-1) variables
%   M12_{i,j,a,b,k}:  (Nx+1)*Nx^3*(Nt-1) variables  (Nx+1 in the 2nd slot)
%   M21_{i,j,a,b,k}:  (Nx+1)*Nx^3*(Nt-1) variables  (Nx+1 in the 3rd slot)
%   M22_{i,j,a,b,k}:  (Nx+1)*Nx^3*(Nt-1) variables  (Nx+1 in the 4th slot)
%

%% Grid (shared by all four coordinates)
x  = linspace(xa,xb,Nx+1)';          % face points, (Nx+1)x1
dx = x(2)-x(1);  dt = 1/(Nt-1);
xhat = 0.5*(x(1:end-1)+x(2:end));    % cell midpoints, Nx x1

%% Marginal densities (evaluated at cell-centre pairs, then normalised)
[XH, YH] = ndgrid(xhat, xhat);                 % Nx x Nx grids
mu1 = arrayfun(@mu1_func, XH, YH);  mu1 = mu1/sum(mu1(:));   % Nx x Nx
mu2 = arrayfun(@mu2_func, XH, YH);  mu2 = mu2/sum(mu2(:));   % Nx x Nx

%% Variable sizes
NGamma = Nx^4*Nt;
NM11   = (Nx+1)*Nx^3*(Nt-1);
NM12   = (Nx+1)*Nx^3*(Nt-1);
NM21   = (Nx+1)*Nx^3*(Nt-1);
NM22   = (Nx+1)*Nx^3*(Nt-1);
Ntot   = NGamma + NM11 + NM12 + NM21 + NM22;

%% Constraints and objective matrix
[Aeq, beq] = constraints_2d(x, mu1, mu2, Nt, Nx);

% D_mat_gen_2d already bakes the (z1-z2)-dependent weights into the
% matrix entries (see that file for why this cannot be a post-hoc
% diagonal rescaling in 2-D, unlike the 1-D code's s_vec/spdiags step).
D_mat = D_mat_gen_2d(x, Nt, Nx);   % (Nx^4*(Nt-1)) x Ntot

%% Bounds: Gamma >= 0; all M components unconstrained
lb           = -Inf*ones(Ntot,1);
lb(1:NGamma) = 0;
ub           = Inf*ones(Ntot,1);

%% Solve via LP reformulation: min sum(t) s.t. D*sol - t <= 0, -D*sol - t <= 0
[sol, fval] = solve_abs_diff_linprog(D_mat, Aeq, beq, lb, ub);

%% Evaluate objective and extract solution arrays
const_error = Aeq*sol - beq;

% Objective: sum(|D_mat*sol|) * dt   (weights already in D_mat; Gamma/M
% are cell MASSES here, matching mu1/mu2's normalization to sum to 1 --
% not densities -- so no extra dx^4 factor is needed, exactly as the
% 1-D code only multiplies by dt and not by dx*dy.)
cost    = zeros(2,1);
cost(1) = sum(abs(D_mat*sol)) * dt;

% Reshape solution
Pi  = reshape(sol(1:NGamma), Nx, Nx, Nx, Nx, Nt);
% (M components left unreshaped here; unpack similarly if needed:
%  M11 = reshape(sol(NGamma+1:NGamma+NM11), Nx+1, Nx, Nx, Nx, Nt-1); etc.)

Pi_final = reshape(Pi(:,:,:,:,end), Nx*Nx, Nx*Nx);   % coupling at t=1,

% Alternative cost: W2^2 from the terminal coupling directly
% dist_sq(i,j,a,b) = (xhat_i-xhat_a)^2 + (xhat_j-xhat_b)^2
[I,J,A,B] = ndgrid(1:Nx,1:Nx,1:Nx,1:Nx);
dist_sq = (xhat(I(:))-xhat(A(:))).^2 + (xhat(J(:))-xhat(B(:))).^2;
Pi_end  = Pi(:,:,:,:,end);
cost(2) = sum( Pi_end(:) .* dist_sq );

%% Diagnostics
mu1_marg = sum(sum(Pi_end, 4), 3);   % Nx x Nx, sum over (a,b)
mu2_marg = squeeze(sum(sum(Pi_end, 1), 2));  % Nx x Nx, sum over (i,j)

diag.const_err_max   = max(abs(const_error));
diag.row_sum_err_max = max(abs(mu1_marg(:) - mu1(:)));
diag.col_sum_err_max = max(abs(mu2_marg(:) - mu2(:)));

end