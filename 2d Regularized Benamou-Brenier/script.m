clear; clc; Nt_list = [2; 4; 8]; N_list = [2; 4; 8; 16];
cost = zeros(4,3); err = zeros(4,3); timeval = zeros(4,3);

cost_ex = 0.013321449188603; xa = 0; xb = 1; % Ex 1
% cost_ex = 0.015511141991097; xa = -0.5; xb = 0.5; % Ex 2

beta = 10^(-3); maxit = 40; alpha = 0.3;
for i=1:4
    for j=1:3
        n = N_list(i)+1; dx = 1/(n-1);
        rho0 = zeros(n); rho1 = zeros(n);
        x1 = linspace(xa,xb,n); x2 = linspace(xa,xb,n);

        for k=1:n
            for l=1:n
                rho0(k,l) = rho0_func(x1(k),x2(l));
                rho1(k,l) = rho1_func(x1(k),x2(l));
            end
        end
        rho0 = rho0/sum(rho0(:)); rho1 = rho1/sum(rho1(:));

        tic; cost(i,j) = quad_solver(rho0, rho1, n, Nt_list(j)-1, dx, beta, maxit, alpha, xa, xb);
        timeval(i,j) = toc; err(i,j) = abs(cost(i,j) - cost_ex);
    end
end