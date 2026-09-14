function [cost, m, k, rho, f_val] = quad_solver(Nx, Nt, beta, maxit, alpha, xa, xb)

clear options;

%% Grid
x  = linspace(xa,xb,Nx+1)';          % z1 face points, (Nx+1)×1
y  = linspace(xa,xb,Nx+1)';          % z2 face points, (Ny+1)×1

xhat = 0.5*(x(1:end-1)+x(2:end));  % z1 cell midpoints, Nx×1
yhat = 0.5*(y(1:end-1)+y(2:end));  % z2 cell midpoints, Ny×1

dx = xhat(2)-xhat(1); % total_matches = sum(ismember(xhat, yhat)); [size(xhat,1) size(yhat,1) total_matches]

%% Marginal densities (evaluated at cell midpoints, then normalised)
mu1 = arrayfun(@rho0_func, xhat); mu1 = mu1+1e-6;  mu1 = mu1/sum(mu1);   % Nx×1
mu2 = arrayfun(@rho1_func, yhat); mu2 = mu2+1e-6;  mu2 = mu2/sum(mu2);   % Ny×1

% Initial Guess
[m, rho] = iter0(mu1, mu2, Nx, Nt, dx);
u = [m(:);rho(:)]; dt = 1/(Nt+1); f_val = zeros(maxit,1);
n1 = (Nx-1)*(Nt+1); n2 = Nx*Nt; N = n1+n2; M = (Nt+1)*(Nx-1)+Nt;

% Defining matrix A and its null space
[A,~] = gen_A(Nx, Nt, mu1, mu2, dx);
[~,~,V] = svd(A); W = V(:,M+1:N); W2 = W(n1+1:end,:);

for k=1:maxit
    
    u0 = u; [m, rho] = vet2mat(Nx, Nt, u0); u_bar_2 = u(n1+1:end);
    f_val(k) = func_val(Nx, Nt, beta, m, rho, mu1, dx);

    f_grad = grad_full(Nx, Nt, beta, m, rho, mu1, dx);
    hess_val = hess_full(Nx, Nt, beta, m, rho, mu1);
    hessian = W'*hess_val*W; hessian = (hessian + hessian')/2;

    c = quadprog(hessian,W'*f_grad,-W2,u_bar_2,[],[],[],[],[]);
    u = u0 + alpha*W*c;

    if norm(c, 'Inf') <= 10^(-12)
        f_val = f_val(1:k,1);
        break;
    end

end

f_cost = @(x) rho1_func(x).*log(rho1_func(x)) - rho0_func(x).*log(rho0_func(x));
cost = dt*func_val(Nx, Nt, beta, m, rho, mu1, dx) + 2*beta*integral(f_cost,xa,xb);

end
