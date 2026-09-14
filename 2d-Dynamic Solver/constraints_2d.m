function [A, b] = constraints_2d(x, mu1, mu2, Nt, Nx)
% CONSTRAINTS_2D  Build equality constraint matrix A and rhs b for the
%                 dynamic OT problem (flux formulation) when z1, z2 are
%                 each points in R^2 (2-D marginals mu1, mu2 on [0,1]^2).
%
% Variable layout in the solution vector [Gamma; M11; M12; M21; M22]:
%
%   Gamma_{i,j,a,b,k}  at cell centres (xhat_i,xhat_j,xhat_a,xhat_b,t_k)
%                      i,j,a,b = 1..Nx,  k = 1..Nt
%                      z1 = (xhat_i,xhat_j),  z2 = (xhat_a,xhat_b)
%
%   M11_{i,j,a,b,k}   flux of z1's 1st coord, at z1^1-faces (x_i,xhat_j,xhat_a,xhat_b,t_k)
%                      i = 1..Nx+1, j,a,b = 1..Nx, k = 1..Nt-1
%   M12_{i,j,a,b,k}   flux of z1's 2nd coord, at z1^2-faces (xhat_i,x_j,xhat_a,xhat_b,t_k)
%                      j = 1..Nx+1, i,a,b = 1..Nx, k = 1..Nt-1
%   M21_{i,j,a,b,k}   flux of z2's 1st coord, at z2^1-faces (xhat_i,xhat_j,x_a,xhat_b,t_k)
%                      a = 1..Nx+1, i,j,b = 1..Nx, k = 1..Nt-1
%   M22_{i,j,a,b,k}   flux of z2's 2nd coord, at z2^2-faces (xhat_i,xhat_j,xhat_a,x_b,t_k)
%                      b = 1..Nx+1, i,j,a = 1..Nx, k = 1..Nt-1
%
% All five variables use MATLAB's column-major linear index over the
% 5-tuple (i,j,a,b,k) via sub2ind -- this generalises the i-fastest,
% then-j, ..., then-k convention of the 1-D code to four spatial indices.
%
% x         - shared 1-D face grid (Nx+1)x1 on [0,1], used for ALL FOUR
%             coordinate directions (both densities live on the unit square).
% mu1, mu2  - Nx x Nx matrices, mu1(i,j) = mass of mu1 in cell (xhat_i,xhat_j),
%             mu2(a,b) = mass of mu2 in cell (xhat_a,xhat_b); both already
%             normalised to sum to 1.
% Nt, Nx    - number of time points / cells per spatial direction.

dt = 1/(Nt-1);
dx = x(2) - x(1);

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

idxG   = @(i,j,a,b,k)         sub2ind(sizeG,   i,j,a,b,k);
idxM11 = @(i,j,a,b,k) offM11 + sub2ind(sizeM11, i,j,a,b,k);
idxM12 = @(i,j,a,b,k) offM12 + sub2ind(sizeM12, i,j,a,b,k);
idxM21 = @(i,j,a,b,k) offM21 + sub2ind(sizeM21, i,j,a,b,k);
idxM22 = @(i,j,a,b,k) offM22 + sub2ind(sizeM22, i,j,a,b,k);

[I,J,A,B] = ndgrid(1:Nx,1:Nx,1:Nx,1:Nx);
Iv = I(:); Jv = J(:); Av = A(:); Bv = B(:);
n_cell = numel(Iv);   % = Nx^4

%% -----------------------------------------------------------------------
%% (D1) Initial condition: Gamma_{i,j,a,b,1} = p_{i,j,a,b}
%       p is the reference measure on the discrete diagonal
%       {(i,j) == (a,b)}: mass 1/Nx^2 on each of the Nx^2 diagonal cells.
%       Nx^4 constraints (one per cell, at k=1).
%% -----------------------------------------------------------------------
n_D1 = n_cell;
rows_D1 = (1:n_D1)';
cols_D1 = idxG(Iv, Jv, Av, Bv, ones(n_D1,1));
vals_D1 = ones(n_D1,1);

diag_mask = (Iv == Av) & (Jv == Bv);
n_diag    = sum(diag_mask);          % = Nx^2
b_D1 = zeros(n_D1,1);
b_D1(diag_mask) = 1/n_diag;

%% -----------------------------------------------------------------------
%% (D2) Continuity equation: Nx^4*(Nt-1) constraints
%
%  (Gamma_{i,j,a,b,k+1}-Gamma_{i,j,a,b,k})/dt
%  + (M11_{i+1,j,a,b,k}-M11_{i,j,a,b,k})/dx
%  + (M12_{i,j+1,a,b,k}-M12_{i,j,a,b,k})/dx
%  + (M21_{i,j,a+1,b,k}-M21_{i,j,a,b,k})/dx
%  + (M22_{i,j,a,b+1,k}-M22_{i,j,a,b,k})/dx = 0
%% -----------------------------------------------------------------------
n_D2  = n_cell*(Nt-1);
r2 = zeros(n_D2*10,1); c2 = zeros(n_D2*10,1); v2 = zeros(n_D2*10,1);
ptr = 0;
row0 = n_D1;

for k = 1:Nt-1
    rows = row0 + (k-1)*n_cell + (1:n_cell)';
    kk  = k*ones(n_cell,1);
    kk1 = (k+1)*ones(n_cell,1);

    cols = [ idxG(Iv,Jv,Av,Bv,kk),   idxG(Iv,Jv,Av,Bv,kk1), ...
             idxM11(Iv,  Jv,Av,Bv,kk), idxM11(Iv+1,Jv,Av,Bv,kk), ...
             idxM12(Iv,Jv,  Av,Bv,kk), idxM12(Iv,Jv+1,Av,Bv,kk), ...
             idxM21(Iv,Jv,Av,  Bv,kk), idxM21(Iv,Jv,Av+1,Bv,kk), ...
             idxM22(Iv,Jv,Av,Bv,  kk), idxM22(Iv,Jv,Av,Bv+1,kk) ];

    vals = repmat([-1/dt, 1/dt, -1/dx, 1/dx, -1/dx, 1/dx, -1/dx, 1/dx, -1/dx, 1/dx], n_cell, 1);
    rows_rep = repmat(rows, 1, 10);

    r2(ptr+1:ptr+10*n_cell) = rows_rep(:);
    c2(ptr+1:ptr+10*n_cell) = cols(:);
    v2(ptr+1:ptr+10*n_cell) = vals(:);
    ptr = ptr + 10*n_cell;
end

%% -----------------------------------------------------------------------
%% (D3) Terminal marginal constraints at k = Nt:  2*Nx^2 constraints
%
%   z1-marginal:  sum_{a,b} Gamma_{i,j,a,b,Nt} = mu1(i,j),  i,j = 1..Nx
%   z2-marginal:  sum_{i,j} Gamma_{i,j,a,b,Nt} = mu2(a,b),  a,b = 1..Nx
%% -----------------------------------------------------------------------
row_offset = n_D1 + n_D2;

row_z1 = row_offset + (Iv + (Jv-1)*Nx);          % Nx^2 distinct rows, repeated
row_z2 = row_offset + Nx^2 + (Av + (Bv-1)*Nx);    % Nx^2 distinct rows, repeated

col_term = idxG(Iv,Jv,Av,Bv, Nt*ones(n_cell,1));

r3 = [row_z1; row_z2];
c3 = [col_term; col_term];
v3 = ones(2*n_cell,1);

b_D3 = [mu1(:); mu2(:)];   % column-major: matches i+(j-1)*Nx / a+(b-1)*Nx ordering

%% -----------------------------------------------------------------------
%% (D4) No-flux boundary conditions: 8*Nx^3*(Nt-1) constraints
%
%   M11_{1,j,a,b,k} = M11_{Nx+1,j,a,b,k} = 0
%   M12_{i,1,a,b,k} = M12_{i,Nx+1,a,b,k} = 0
%   M21_{i,j,1,b,k} = M21_{i,j,Nx+1,b,k} = 0
%   M22_{i,j,a,1,k} = M22_{i,j,a,Nx+1,k} = 0
%% -----------------------------------------------------------------------
[P,Q,R] = ndgrid(1:Nx,1:Nx,1:Nx);
Pv = P(:); Qv = Q(:); Rv = R(:);
n_face = numel(Pv);   % = Nx^3

n_D3 = 2*Nx^2;
n_D4 = 8*n_face*(Nt-1);
r4 = zeros(n_D4,1); c4 = zeros(n_D4,1);
ptr4 = 0;
row_offset4 = n_D1 + n_D2 + n_D3;

for k = 1:Nt-1
    kk = k*ones(n_face,1);

    % M11 boundary: fix i = 1 or i = Nx+1; free (j,a,b) = (Pv,Qv,Rv)
    blk = row_offset4 + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM11(ones(n_face,1),   Pv,Qv,Rv,kk); ptr4 = ptr4+n_face;
    blk = row_offset4 + n_face*(Nt-1) + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM11((Nx+1)*ones(n_face,1),Pv,Qv,Rv,kk); ptr4 = ptr4+n_face;
end
off = row_offset4 + 2*n_face*(Nt-1);
for k = 1:Nt-1
    kk = k*ones(n_face,1);
    % M12 boundary: fix j = 1 or j = Nx+1; free (i,a,b) = (Pv,Qv,Rv)
    blk = off + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM12(Pv, ones(n_face,1),   Qv,Rv,kk); ptr4 = ptr4+n_face;
    blk = off + n_face*(Nt-1) + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM12(Pv, (Nx+1)*ones(n_face,1),Qv,Rv,kk); ptr4 = ptr4+n_face;
end
off = off + 2*n_face*(Nt-1);
for k = 1:Nt-1
    kk = k*ones(n_face,1);
    % M21 boundary: fix a = 1 or a = Nx+1; free (i,j,b) = (Pv,Qv,Rv)
    blk = off + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM21(Pv,Qv, ones(n_face,1),   Rv,kk); ptr4 = ptr4+n_face;
    blk = off + n_face*(Nt-1) + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM21(Pv,Qv, (Nx+1)*ones(n_face,1),Rv,kk); ptr4 = ptr4+n_face;
end
off = off + 2*n_face*(Nt-1);
for k = 1:Nt-1
    kk = k*ones(n_face,1);
    % M22 boundary: fix b = 1 or b = Nx+1; free (i,j,a) = (Pv,Qv,Rv)
    blk = off + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM22(Pv,Qv,Rv, ones(n_face,1)   ,kk); ptr4 = ptr4+n_face;
    blk = off + n_face*(Nt-1) + (k-1)*n_face + (1:n_face)';
    r4(ptr4+1:ptr4+n_face) = blk; c4(ptr4+1:ptr4+n_face) = idxM22(Pv,Qv,Rv, (Nx+1)*ones(n_face,1),kk); ptr4 = ptr4+n_face;
end
v4 = ones(n_D4,1);

%% -----------------------------------------------------------------------
%% Assemble sparse A and rhs b
%% -----------------------------------------------------------------------
n_rows = n_D1 + n_D2 + n_D3 + n_D4;

rows_all = [rows_D1; r2; r3; r4];
cols_all = [cols_D1; c2; c3; c4];
vals_all = [vals_D1; v2; v3; v4];

A = sparse(rows_all, cols_all, vals_all, n_rows, Ntot);

b = [b_D1; zeros(n_D2,1); b_D3; zeros(n_D4,1)];

end