function [A, b] = constraints(x, y, mu1, mu2, Nt, Nx, Ny)
% CONSTRAINTS  Build equality constraint matrix A and rhs b for the
%              dynamic OT problem (flux formulation, quadratic cost).
%
% Variable layout in the solution vector [Gamma; M1; M2]:
%
%   Gamma_{i,j,k}  at cell centres (xhat_i, yhat_j, t_k)
%                  i=1..Nx,   j=1..Ny,   k=1..Nt
%                  size: NGamma = Nx*Ny*Nt
%                  linear index: i + (j-1)*Nx + (k-1)*Nx*Ny
%
%   M1_{i,j,k}    at z1-faces (x_i, yhat_j, t_k)
%                  i=1..Nx+1, j=1..Ny,   k=1..Nt-1
%                  size: NM1 = (Nx+1)*Ny*(Nt-1)
%                  linear index: i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*Ny
%
%   M2_{i,j,k}    at z2-faces (xhat_i, y_j, t_k)
%                  i=1..Nx,   j=1..Ny+1, k=1..Nt-1
%                  size: NM2 = Nx*(Ny+1)*(Nt-1)
%                  linear index: i + (j-1)*Nx + (k-1)*Nx*(Ny+1)
%
% All indexing follows x->y->t order (x varies fastest, MATLAB column-major).
% Notation matches LaTeX: Nx cells in z1, Ny cells in z2, Nt time points.

dt = 1/(Nt-1);  dx = x(2)-x(1);  dy = y(2)-y(1);
xhat = 0.5*(x(1:end-1)+x(2:end));  % z1 cell midpoints, Nx×1
yhat = 0.5*(y(1:end-1)+y(2:end));  % z2 cell midpoints, Ny×1

NGamma = Nx*Ny*Nt;
NM1    = (Nx+1)*Ny*(Nt-1);
NM2    = Nx*(Ny+1)*(Nt-1);
Ntot   = NGamma + NM1 + NM2;

% Helper: linear index functions (return column vector of indices)
idxG  = @(i,j,k)  i + (j-1)*Nx     + (k-1)*Nx*Ny;          % Gamma
idxM1 = @(i,j,k)  i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*Ny;      % M1
idxM2 = @(i,j,k)  i + (j-1)*Nx     + (k-1)*Nx*(Ny+1);      % M2

%% -----------------------------------------------------------------------
%% (D1)  Initial condition:  Gamma_{i,j,1} = p_{i,j}
%         p is the diagonal reference measure: 1/min(Nx,Ny) if i==j, else 0.
%         This gives Nx*Ny equality constraints (one per cell at k=1).
% -----------------------------------------------------------------------
n_D1 = Nx*Ny;
rows_D1 = (1:n_D1)';
cols_D1 = idxG( reshape(repmat(1:Nx,Ny,1),[],1), ...
                repmat((1:Ny)',Nx,1), ones(n_D1,1) );
vals_D1 = ones(n_D1,1);

b_D1 = zeros(n_D1,1);
for i = 1:min(Nx,Ny)
    b_D1( idxG(i,i,1) ) = 1/min(Nx,Ny);
end

% b_D1 = zeros(n_D1, 1); tol = 1e-12;
% for i = 1:Nx
%     for j = 1:Ny
%         if abs(xhat(i) - yhat(j)) < tol
%             b_D1(idxG(i, j, 1)) = 1;
%         end
%     end
% end
% b_D1 = b_D1 / sum(b_D1);

%% (D1) Initial condition: p on physical diagonal {z1 = z2}
% b_D1 = zeros(n_D1, 1);
% diag_cells = [];
% for i = 1:Nx
%     [dist, j_nearest] = min(abs(xhat(i) - yhat));
%     if dist < 0.5*max(dx, dy)
%         diag_cells(end+1,:) = [i, j_nearest];
%     end
% end
% N_diag = size(diag_cells, 1);
% for k = 1:N_diag
%     b_D1(idxG(diag_cells(k,1), diag_cells(k,2), 1)) = 1/N_diag;
% end

%% -----------------------------------------------------------------------
%% (D2)  Continuity equation (Nt-1)*Nx*Ny constraints:
%
%   (Gamma_{i,j,k+1} - Gamma_{i,j,k})/dt
%   + (M1_{i+1,j,k}  - M1_{i,j,k}  )/dx
%   + (M2_{i,j+1,k}  - M2_{i,j,k}  )/dy  = 0
%
%   for i=1..Nx, j=1..Ny, k=1..Nt-1
% -----------------------------------------------------------------------
n_D2   = Nx*Ny*(Nt-1);
n_nnz  = n_D2*6;               % 2 Gamma + 2 M1 + 2 M2 per row
r2 = zeros(n_nnz,1);  c2 = zeros(n_nnz,1);  v2 = zeros(n_nnz,1);
ptr = 0;

for k = 1:Nt-1
    for j = 1:Ny
        for i = 1:Nx
            row = n_D1 + idxG(i,j,k);   % row offset within D2 block

            % Gamma time derivative
            entries = [ row, idxG(i,j,k),   -1/dt;
                        row, idxG(i,j,k+1),  1/dt;
                        row, NGamma + idxM1(i,  j,k), -1/dx;
                        row, NGamma + idxM1(i+1,j,k),  1/dx;
                        row, NGamma + NM1 + idxM2(i,j,  k), -1/dy;
                        row, NGamma + NM1 + idxM2(i,j+1,k),  1/dy ];

            r2(ptr+1:ptr+6) = entries(:,1);
            c2(ptr+1:ptr+6) = entries(:,2);
            v2(ptr+1:ptr+6) = entries(:,3);
            ptr = ptr + 6;
        end
    end
end

%% -----------------------------------------------------------------------
%% (D3)  Terminal marginal constraints at k=Nt:
%
%   z1-marginal:  sum_{j=1}^{Ny} Gamma_{i,j,Nt} = mu1(xhat_i),  i=1..Nx
%   z2-marginal:  sum_{i=1}^{Nx} Gamma_{i,j,Nt} = mu2(yhat_j),  j=1..Ny
% -----------------------------------------------------------------------
n_D3 = Nx + Ny;
r3 = zeros(n_D3*max(Nx,Ny),1);
c3 = zeros(n_D3*max(Nx,Ny),1);
v3 = zeros(n_D3*max(Nx,Ny),1);
ptr3 = 0;
row_offset = n_D1 + n_D2;

% z1-marginal: row per i, sum over j
for i = 1:Nx
    for j = 1:Ny
        ptr3 = ptr3+1;
        r3(ptr3) = row_offset + i;
        c3(ptr3) = idxG(i,j,Nt);
        v3(ptr3) = 1;
    end
end
% z2-marginal: row per j, sum over i
for j = 1:Ny
    for i = 1:Nx
        ptr3 = ptr3+1;
        r3(ptr3) = row_offset + Nx + j;
        c3(ptr3) = idxG(i,j,Nt);
        v3(ptr3) = 1;
    end
end
r3 = r3(1:ptr3);  c3 = c3(1:ptr3);  v3 = v3(1:ptr3);

b_D3 = [mu1; mu2];    % mu1 is Nx×1, mu2 is Ny×1

%% -----------------------------------------------------------------------
%% (D4)  No-flux boundary conditions:
%
%   M1_{1,  j,k} = 0   and   M1_{Nx+1,j,k} = 0,   j=1..Ny,  k=1..Nt-1
%   M2_{i,1,  k} = 0   and   M2_{i,Ny+1,k} = 0,   i=1..Nx,  k=1..Nt-1
% -----------------------------------------------------------------------
n_D4_M1 = 2*Ny*(Nt-1);
n_D4_M2 = 2*Nx*(Nt-1);
n_D4    = n_D4_M1 + n_D4_M2;
r4 = zeros(n_D4,1);  c4 = zeros(n_D4,1);  v4 = ones(n_D4,1);
ptr4 = 0;
row_offset4 = n_D1 + n_D2 + n_D3;

for k = 1:Nt-1
    for j = 1:Ny
        % M1 left boundary  (i=1)
        ptr4 = ptr4+1;
        r4(ptr4) = row_offset4 + (k-1)*Ny + j;
        c4(ptr4) = NGamma + idxM1(1,j,k);
        % M1 right boundary (i=Nx+1)
        ptr4 = ptr4+1;
        r4(ptr4) = row_offset4 + Ny*(Nt-1) + (k-1)*Ny + j;
        c4(ptr4) = NGamma + idxM1(Nx+1,j,k);
    end
end
row_offset4b = row_offset4 + n_D4_M1;
for k = 1:Nt-1
    for i = 1:Nx
        % M2 bottom boundary (j=1)
        ptr4 = ptr4+1;
        r4(ptr4) = row_offset4b + (k-1)*Nx + i;
        c4(ptr4) = NGamma + NM1 + idxM2(i,1,k);
        % M2 top boundary (j=Ny+1)
        ptr4 = ptr4+1;
        r4(ptr4) = row_offset4b + Nx*(Nt-1) + (k-1)*Nx + i;
        c4(ptr4) = NGamma + NM1 + idxM2(i,Ny+1,k);
    end
end

%% -----------------------------------------------------------------------
%% Assemble sparse A and rhs b
% -----------------------------------------------------------------------
n_rows = n_D1 + n_D2 + n_D3 + n_D4;

rows_all = [rows_D1;          r2(1:ptr);  r3;  r4];
cols_all = [cols_D1;          c2(1:ptr);  c3;  c4];
vals_all = [vals_D1; v2(1:ptr); v3; ones(length(r4),1)];

A = sparse(rows_all, cols_all, vals_all, n_rows, Ntot);

b = [b_D1;
     zeros(n_D2, 1);
     b_D3;
     zeros(n_D4, 1)];

end