function [A,b] = gen_A(n, L, rho0, rho1, dx)

%(L+1)*(n-1)+L - number of equations; (n-1)*(L+1) + n*L - number of unknowns
A = zeros((L+1)*(n-1)+L,(n-1)*(L+1) + n*L);
dt = 1/(L+1); b = zeros((L+1)*(n-1)+L,1);

% Defining matrix A (Forward Difference Used)
for l=0:L
    for i = 1:n-1
        k = i + (n-1)*l;
        % rho coef
        if l==0
            b(k) = rho0(i)/dt;
            A(k, (L+1)*(n-1) + l*n + i) = 1/dt;
        elseif l==L
            b(k) = -rho1(i)/dt;
            A(k, (L+1)*(n-1) + (l-1)*n + i) = -1/dt;
        else
            A(k, (L+1)*(n-1) + (l-1)*n + i) = -1/dt;
            A(k, (L+1)*(n-1) + l*n + i) = 1/dt;
        end
        % m coef
        i_x = i-1;
        if i==1
            A(k, l*(n-1)+i_x+1) = 1/dx;
        else
            A(k, l*(n-1)+i_x) = -1/dx;
            A(k, l*(n-1)+i_x+1) = 1/dx;
        end
    end
end

for l=1:L
    k = (L+1)*(n-1)+l;
    for i=1:n
        A(k, (L+1)*(n-1) + (l-1)*n + i) = 1;
    end
    b(k) = 1;
end

end