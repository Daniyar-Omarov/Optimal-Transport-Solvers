function D_mat = D_mat_gen(Nt, Nx, Ny)
% D_MAT_GEN  Build the difference matrix D for the objective function.
%
% The objective (flux formulation, corrected) is:
%
%   sum_{k=1}^{Nt-1} sum_{i=1}^{Nx} sum_{j=1}^{Ny}
%       s_{i,j} * | avg_M1_{i,j,k} - avg_M2_{i,j,k} | * dx*dy*dt
%
% where
%   avg_M1_{i,j,k} = (M1_{i,j,k} + M1_{i+1,j,k}) / 2   (M1 interpolated to cell centre)
%   avg_M2_{i,j,k} = (M2_{i,j,k} + M2_{i,j+1,k}) / 2   (M2 interpolated to cell centre)
%   s_{i,j}        = 2 * |xhat_i - yhat_j|               (weight from nabla c)
%
% D_mat encodes the linear map sol -> (avg_M1 - avg_M2) evaluated at every
% (i,j,k) triple, so that the objective equals
%   s .* |D_mat * sol| summed with measure weights dx*dy*dt
% (the s weights and dx*dy*dt scaling are applied in lin_min, not here).
%
% Variable layout in sol = [Gamma; M1; M2]  (x->y->t order):
%
%   Gamma_{i,j,k}: index i + (j-1)*Nx     + (k-1)*Nx*Ny        (i=1..Nx, j=1..Ny, k=1..Nt)
%   M1_{i,j,k}:   index i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*Ny    (i=1..Nx+1, j=1..Ny, k=1..Nt-1)
%   M2_{i,j,k}:   index i + (j-1)*Nx     + (k-1)*Nx*(Ny+1)    (i=1..Nx, j=1..Ny+1, k=1..Nt-1)

NGamma = Nx*Ny*Nt;
NM1    = (Nx+1)*Ny*(Nt-1);
% NM2  = Nx*(Ny+1)*(Nt-1);   % not needed here but listed for clarity

idxM1 = @(i,j,k)  i + (j-1)*(Nx+1) + (k-1)*(Nx+1)*Ny;
idxM2 = @(i,j,k)  i + (j-1)*Nx     + (k-1)*Nx*(Ny+1);

n_rows = Nx*Ny*(Nt-1);    % one row per (i,j,k), k=1..Nt-1
n_nnz  = n_rows*4;        % 2 M1 entries + 2 M2 entries per row

rows = zeros(n_nnz,1);
cols = zeros(n_nnz,1);
vals = zeros(n_nnz,1);
ptr  = 0;

for k = 1:Nt-1
    for j = 1:Ny
        for i = 1:Nx
            row = i + (j-1)*Nx + (k-1)*Nx*Ny;   % row index in D_mat

            % +1/2 * M1_{i,  j,k}
            ptr = ptr+1;
            rows(ptr) = row;
            cols(ptr) = NGamma + idxM1(i,  j,k);
            vals(ptr) = 1/2;

            % +1/2 * M1_{i+1,j,k}
            ptr = ptr+1;
            rows(ptr) = row;
            cols(ptr) = NGamma + idxM1(i+1,j,k);
            vals(ptr) = 1/2;

            % -1/2 * M2_{i,j,  k}
            ptr = ptr+1;
            rows(ptr) = row;
            cols(ptr) = NGamma + NM1 + idxM2(i,j,  k);
            vals(ptr) = -1/2;

            % -1/2 * M2_{i,j+1,k}
            ptr = ptr+1;
            rows(ptr) = row;
            cols(ptr) = NGamma + NM1 + idxM2(i,j+1,k);
            vals(ptr) = -1/2;
        end
    end
end

Ntot  = NGamma + NM1 + Nx*(Ny+1)*(Nt-1);
D_mat = sparse(rows, cols, vals, n_rows, Ntot);

end