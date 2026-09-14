function value = hess_full(n, L, dx, beta, m_x, m_y, rho, rho0)

value1 = zeros(2*n*(n-1)*(L+1),1);
for l=0:L
    if l==0
        for i=1:n-1
            for j=1:n
                value1((j-1)*(n-1)+i) = 4/(rho0(i,j)+rho0(i+1,j));
            end
        end
        for i=1:n
            for j=1:n-1
                value1((L+1)*n*(n-1)+(j-1)*n+i) = 4/(rho0(i,j)+rho0(i,j+1));
            end
        end
    else
        for i=1:n-1
            for j=1:n
                value1(l*(n-1)*n+(j-1)*(n-1)+i) = 4/(rho(i,j,l)+rho(i+1,j,l));
            end
        end
        for i=1:n
            for j=1:n-1
                value1((L+1)*n*(n-1)+l*(n-1)*n+(j-1)*n+i) = 4/(rho(i,j,l)+rho(i,j+1,l));
            end
        end
    end
    
end

value2 = zeros((n-1)*n*(L+1),L*n^2);
for l = 1:L
    for j=1:n
        for i=1:n-1
            value2((j-1)*(n-1)+i+l*n*(n-1), (j-1)*n+i+(l-1)*n*n) = -4*m_x(i,j,l+1)/(rho(i,j,l)+rho(i+1,j,l))^2;
            value2((j-1)*(n-1)+i+l*n*(n-1), (j-1)*n+i+1+(l-1)*n*n) = -4*m_x(i,j,l+1)/(rho(i,j,l)+rho(i+1,j,l))^2;
        end
    end
end

value3 = zeros((n-1)*n*(L+1),L*n^2);
for l = 1:L
    for j=1:n-1
        for i=1:n
            value3((j-1)*n+i+l*n*(n-1), (j-1)*n+i+(l-1)*n*n) = -4*m_y(i,j,l+1)/(rho(i,j,l)+rho(i,j+1,l))^2;
            value3((j-1)*n+i+l*n*(n-1), j*n+i+(l-1)*n*n) = -4*m_y(i,j,l+1)/(rho(i,j,l)+rho(i,j+1,l))^2;
        end  
    end
end

value4 = zeros(L*n^2);
for l = 1:L
    value4(1+(l-1)*n^2:l*n^2,1+(l-1)*n^2:l*n^2) = J_rr(rho(:,:,l), m_x(:,:,l+1), m_y(:,:,l+1), n) + (beta^2)*I_rr(rho(:,:,l), n, dx);
end

value = sparse(2*n*(n-1)*(L+1) + L*n^2, 2*n*(n-1)*(L+1) + L*n^2);
value(1:n*(n-1)*(L+1),2*n*(n-1)*(L+1)+1:end) = value2;
value(n*(n-1)*(L+1)+1:2*n*(n-1)*(L+1),2*n*(n-1)*(L+1)+1:end) = value3;
value = value+value';
for i=1:2*n*(n-1)*(L+1)
    value(i,i) = value1(i);
end
value(2*n*(n-1)*(L+1)+1:end,2*n*(n-1)*(L+1)+1:end) = value4;

end