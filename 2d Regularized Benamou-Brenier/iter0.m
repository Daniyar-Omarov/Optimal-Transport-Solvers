function [m_x, m_y, rho] = iter0(rho0, rho1, n, L, dx)

% n - grid size; L+1 - number of interval in time;
% dx - mesh size in x and y; dt - mesh size in time;
dt = 1/(L+1); t_l = dt*linspace(0, L+1, L+2);

rho = zeros(n,n,L);

for l=1:L
    for j=1:n
        for i=1:n
            rho(i,j,l) = (1-t_l(l+1))*rho0(i,j) + t_l(l+1)*rho1(i,j);
        end
    end
end

m_x = zeros(n-1,n,L+1); m_y = zeros(n,n-1,L+1);
for l=0:L
    temp_mat = zeros(n^2, 2*(n-1)*n); temp_vec = zeros(n^2,1);
    for j=1:n
        for i=1:n
            k = i+(j-1)*n;
            if i==1 && j==1
                temp_mat(k, 1) = 1; temp_mat(k,(n-1)*n+1) = 1;
            elseif i==1 && j==n
                temp_mat(k, 1+ (n-1)^2) = 1; temp_mat(k,(n-1)*n+n*(n-2)+1) = -1;
            elseif i==n && j==1
                temp_mat(k, n-1) = -1; temp_mat(k,(n-1)*n+n) = 1;
            elseif i==n && j==n
                temp_mat(k, (n-1)*n) = -1; temp_mat(k,2*(n-1)*n) = -1;
            elseif j==1
                temp_mat(k, i) = 1; temp_mat(k, i-1) = -1; temp_mat(k, (n-1)*n + i) = 1; 
            elseif j==n
                temp_mat(k, i+(n-1)^2) = 1; temp_mat(k, i-1+(n-1)^2) = -1; temp_mat(k, (n-1)*n + i + (n-2)*n) = -1; 
            elseif i==1
                temp_mat(k, 1+(j-1)*(n-1)) = 1; temp_mat(k, (n-1)*n + 1 + (j-1)*n) = 1; temp_mat(k, (n-1)*n + 1 + (j-2)*n) = -1; 
            elseif i==n
                temp_mat(k, n-1 + (j-1)*(n-1)) = -1; temp_mat(k, (n-1)*n + n + (j-1)*n) = 1; temp_mat(k, (n-1)*n + n + (j-2)*n) = -1; 
            else
               temp_mat(k, i+(j-1)*(n-1)) = 1; temp_mat(k, i-1 + (j-1)*(n-1)) = -1; temp_mat(k, (n-1)*n + i + (j-1)*n) = 1; temp_mat(k, (n-1)*n + i + (j-2)*n) = -1; 
            end
            if l==0
                temp_vec(k) = -(dx/dt)*(rho(i,j,l+1)-rho0(i,j));
            elseif l==L
                temp_vec(k) = -(dx/dt)*(rho1(i,j)-rho(i,j,l));
            else
                temp_vec(k) = -(dx/dt)*(rho(i,j,l+1)-rho(i,j,l));
            end
        end
    end

    temp = pinv(temp_mat)*temp_vec;
    m_x(1+l*(n-1)*n:(l+1)*(n-1)*n) = temp(1:n*(n-1));
    m_y(1+l*(n-1)*n:(l+1)*(n-1)*n) = temp(1+n*(n-1):end);
end

end