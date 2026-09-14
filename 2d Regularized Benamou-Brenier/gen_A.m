function [A,b] = gen_A(n, L, dx, rho0, rho1)

dt = 1/(L+1);

% Defining matrix A
% (L+1)*n*n+L - number of equations; 2*n*(n-1)*(L+1) + L*n^2 - number of unknowns
A = zeros((L+1)*n*n+L,2*n*(n-1)*(L+1) + L*n^2);
b = zeros((L+1)*n*n+L,1);
for l=0:L
    for j=1:n
        for i = 1:n
            k = i + (j-1)*n + (n^2)*l;
            % rho coef
            if l==0
                b(k) = dx*rho0(i,j)/dt;
                A(k, 2*n*(n-1)*(L+1) + i + (j-1)*n + (n^2)*l) = dx/dt;
            elseif l==L
                b(k) = -dx*rho1(i,j)/dt;
                A(k, 2*n*(n-1)*(L+1) + i + (j-1)*n + (n^2)*(l-1)) = -dx/dt;
            else
                A(k, 2*n*(n-1)*(L+1) + i + (j-1)*n + (n^2)*l) = dx/dt;
                A(k, 2*n*(n-1)*(L+1) + i + (j-1)*n + (n^2)*(l-1)) = -dx/dt;
            end
            
            % m coef
            if i==1 && j==1
                A(k, 1 + l*n*(n-1)) = 1; A(k,n*(n-1)*(L+1)+1 + l*n*(n-1)) = 1;
            elseif i==1 && j==n
                A(k, 1+ (n-1)^2 + l*n*(n-1)) = 1; A(k,n*(n-1)*(L+1)+n*(n-2)+1 + l*n*(n-1)) = -1;
            elseif i==n && j==1
                A(k, n-1 + l*n*(n-1)) = -1; A(k,n*(n-1)*(L+1)+n + l*n*(n-1)) = 1;
            elseif i==n && j==n
                A(k, (n-1)*n + l*n*(n-1)) = -1; A(k,n*(n-1)*(L+1)+(n-1)*n + l*n*(n-1)) = -1;
            elseif j==1
                A(k, i + l*n*(n-1)) = 1; A(k, i-1 + l*n*(n-1)) = -1; A(k, n*(n-1)*(L+1) + i + l*n*(n-1)) = 1; 
            elseif j==n
                A(k, i+(n-1)^2 + l*n*(n-1)) = 1; A(k, i-1+(n-1)^2 + l*n*(n-1)) = -1; A(k, n*(n-1)*(L+1) + i + (n-2)*n + l*n*(n-1)) = -1; 
            elseif i==1
                A(k, 1+(j-1)*(n-1) + l*n*(n-1)) = 1; A(k, n*(n-1)*(L+1) + 1 + (j-1)*n + l*n*(n-1)) = 1; A(k, n*(n-1)*(L+1) + 1 + (j-2)*n + l*n*(n-1)) = -1; 
            elseif i==n
                A(k, n-1 + (j-1)*(n-1) + l*n*(n-1)) = -1; A(k, n*(n-1)*(L+1) + n + (j-1)*n + l*n*(n-1)) = 1; A(k, n*(n-1)*(L+1) + n + (j-2)*n + l*n*(n-1)) = -1; 
            else
               A(k, i+(j-1)*(n-1) + l*n*(n-1)) = 1; A(k, i-1 + (j-1)*(n-1) + l*n*(n-1)) = -1; A(k, n*(n-1)*(L+1) + i + (j-1)*n + l*n*(n-1)) = 1; A(k, n*(n-1)*(L+1) + i + (j-2)*n + l*n*(n-1)) = -1; 
            end
        end
    end
end

for l=1:L
    k = (L+1)*n*n+l;
    for j=1:n
        for i=1:n
            A(k, 2*n*(n-1)*(L+1) + (l-1)*n^2 + i + (j-1)*n) = 1;
        end
    end
    b(k) = 1;
end

end