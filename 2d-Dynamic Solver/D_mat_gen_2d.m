function D_mat = D_mat_gen_2d(x, Nt, Nx)
% D_MAT_GEN_2D  Build the (already position-weighted) difference matrix D
%               for the 2-D objective function.
%
% The objective (flux formulation, corrected dot-product form) is:
%
%   sum_{k=1}^{Nt-1} sum_{i,j,a,b=1}^{Nx}
%       | (z1-z2) . (avgM1_{i,j,a,b,k} - avgM2_{i,j,a,b,k}) | * (dx)^4 * dt
%
% where z1-z2 = (xhat_i - xhat_a, xhat_j - xhat_b) in R^2,
%       avgM1 = (avg of M11 to centre, avg of M12 to centre),
%       avgM2 = (avg of M21 to centre, avg of M22 to centre),
% and "." is the R^2 Euclidean dot product -- NOT a product of two
% separate norms |z1-z2| * |avgM1-avgM2| (see discussion: those two
% expressions differ whenever avgM1-avgM2 is not parallel to z1-z2).
%
% Expanding the dot product and folding in the factor of 2 from the
% original |nabla c . v| = 2|(z1-z2).(v^(1)-v^(2))| gives, for each row
% (i,j,a,b,k), coefficients (note the 2 * 0.5 = 1 cancellation, so no
% explicit 1/2 appears on the individual face values):
%
%   +dx1 on M11_{i,j,a,b,k}      +dx1 on M11_{i+1,j,a,b,k}
%   -dx1 on M21_{i,j,a,b,k}      -dx1 on M21_{i,j,a+1,b,k}
%   +dx2 on M12_{i,j,a,b,k}      +dx2 on M12_{i,j+1,a,b,k}
%   -dx2 on M22_{i,j,a,b,k}      -dx2 on M22_{i,j,a,b+1,k}
%
% where dx1 = xhat_i - xhat_a,  dx2 = xhat_j - xhat_b.
%
% D_mat encodes exactly this linear map sol -> (row value at every
% (i,j,a,b,k)), so that the objective equals
%   sum( |D_mat * sol| ) * (dx)^4 * dt
% (the (dx)^4 * dt scaling is applied by the caller, e.g. lin_min_2d,
% not here -- matching the convention of the 1-D D_mat_gen.m).
%
% NOTE: unlike the 1-D code, the position weight here cannot be applied
% as a single post-multiplying diagonal (spdiags(s_vec,...)*D_mat),
% because within one row the two position weights dx1, dx2 multiply two
% *different* blocks of columns (the M11/M21 pair vs. the M12/M22 pair).
% So the weights are baked directly into the matrix entries below.

xhat = 0.5*(x(1:end-1)+x(2:end));   % Nx x 1 cell midpoints, shared by all 4 dirs

sizeG   = [Nx,   Nx,   Nx,   Nx,   Nt];
sizeM11 = [Nx+1, Nx,   Nx,   Nx,   Nt-1];
sizeM12 = [Nx,   Nx+1, Nx,   Nx,   Nt-1];
sizeM21 = [Nx,   Nx,   Nx+1, Nx,   Nt-1];
sizeM22 = [Nx,   Nx,   Nx,   Nx+1, Nt-1];

NGamma = prod(sizeG);
NM11   = prod(sizeM11);
NM12   = prod(sizeM12);
NM21   = prod(sizeM21);
NM22   = prod(sizeM22);
Ntot   = NGamma + NM11 + NM12 + NM21 + NM22;

offM11 = NGamma;
offM12 = offM11 + NM11;
offM21 = offM12 + NM12;
offM22 = offM21 + NM21;

idxM11 = @(i,j,a,b,k) offM11 + sub2ind(sizeM11, i,j,a,b,k);
idxM12 = @(i,j,a,b,k) offM12 + sub2ind(sizeM12, i,j,a,b,k);
idxM21 = @(i,j,a,b,k) offM21 + sub2ind(sizeM21, i,j,a,b,k);
idxM22 = @(i,j,a,b,k) offM22 + sub2ind(sizeM22, i,j,a,b,k);

[I,J,A,B] = ndgrid(1:Nx,1:Nx,1:Nx,1:Nx);
Iv = I(:); Jv = J(:); Av = A(:); Bv = B(:);
n_cell = numel(Iv);   % = Nx^4

dx1 = xhat(Iv) - xhat(Av);   % n_cell x 1
dx2 = xhat(Jv) - xhat(Bv);   % n_cell x 1

n_rows = n_cell*(Nt-1);      % one row per (i,j,a,b,k), k = 1..Nt-1
rows = zeros(n_rows*8,1); cols = zeros(n_rows*8,1); vals = zeros(n_rows*8,1);
ptr = 0;

for k = 1:Nt-1
    row_block = (k-1)*n_cell + (1:n_cell)';
    kk = k*ones(n_cell,1);

    cols_k = [ idxM11(Iv,  Jv,Av,Bv,kk), idxM11(Iv+1,Jv,Av,Bv,kk), ...
               idxM21(Iv,Jv,Av,  Bv,kk), idxM21(Iv,Jv,Av+1,Bv,kk), ...
               idxM12(Iv,Jv,  Av,Bv,kk), idxM12(Iv,Jv+1,Av,Bv,kk), ...
               idxM22(Iv,Jv,Av,Bv,  kk), idxM22(Iv,Jv,Av,Bv+1,kk) ];

    vals_k = [ dx1, dx1, -dx1, -dx1, dx2, dx2, -dx2, -dx2 ];

    rows_k = repmat(row_block, 1, 8);

    rows(ptr+1:ptr+8*n_cell) = rows_k(:);
    cols(ptr+1:ptr+8*n_cell) = cols_k(:);
    vals(ptr+1:ptr+8*n_cell) = vals_k(:);
    ptr = ptr + 8*n_cell;
end

D_mat = sparse(rows, cols, vals, n_rows, Ntot);

end