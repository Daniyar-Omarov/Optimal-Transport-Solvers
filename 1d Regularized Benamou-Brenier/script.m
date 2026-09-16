clear; clc; Nt_list = [2; 4; 8]; N_list = [2; 4; 8; 16; 32; 64; 128];
cost = zeros(7,3); timeval = zeros(7,3); err_cost = zeros(7,3);

cost_ex = 1/120; xa = 0; xb = 1; ya = 0; yb = 1; % Ex 1
% cost_ex = (25/16)*asin(4/5)-17/12; xa = 0; xb = 1; ya = 0; yb = 1; % Ex 2
% cost_ex = 1/30; xa = -1; xb = 1; ya = -1; yb = 1; % Ex 3
% cost_ex = 1.25; xa = -10; xb = 10; ya = -10; yb = 10; % Ex 4

beta = 1e-5; alpha = 0.5; maxit = 100;
for i=1:7
    for j=1:3
        tic; cost(i,j) = quad_solver(N_list(i)+1, Nt_list(j), beta, maxit, alpha, xa, xb);
        timeval(i,j) = toc; err_cost(i,j) = abs(cost(i,j) - cost_ex);
    end
end