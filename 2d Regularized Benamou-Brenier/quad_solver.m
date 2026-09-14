function [cost, m_x, m_y, k, rho, f_val] = quad_solver(rho0, rho1, n, L, dx, beta, maxit, alpha, xa, xb)

clear options;
% Initial Guess
[m_x, m_y, rho] = iter0(rho0, rho1, n, L, dx);
u = [m_x(:);m_y(:);rho(:)]; dt = 1/(L+1); f_val = zeros(maxit,1);
n1 = 2*n*(n-1)*(L+1); n2 = L*n^2; N = n1+n2; M = (L+1)*n*n+L;

% Defining matrix A and its null space
[A,~] = gen_A(n, L, dx, rho0, rho1);
[Q,~] = qr(A'); V1 = Q(:,M-(L+1)+1:N); W = V1./vecnorm(V1); W2 = W(n1+1:end,:);

for k=1:maxit
    u0 = u; [m_x, m_y, rho] = vet2mat(n, L, u0); u_bar_2 = u(n1+1:end);

    f_val(k) = func_val(n, L, dx, beta, m_x, m_y, rho, rho0);
    f_grad = grad_full(n, L, dx, beta, m_x, m_y, rho, rho0);
    hess_val = hess_full(n, L, dx, beta, m_x, m_y, rho, rho0);
    hessian = full(sparse(W')*sparse(hess_val)*sparse(W)); hessian = (hessian + hessian')/2;

    options = mskoptimset('Display','off');
    c = quadprog(hessian,W'*f_grad,-W2,u_bar_2,[],[],[],[],[], options);

    u = u0 + alpha*W*c;
    if norm(c, 'Inf') <= 10^(-6)
        f_val = f_val(1:k,1);
        break;
    end

end

cost = dt*func_val(n, L, dx, beta, m_x, m_y, rho, rho0) + 2*beta*simpson2d(n, xa, xb);

end
