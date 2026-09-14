function [m_x, m_y, rho] = vet2mat(n, L, u)

m_x = zeros(n-1,n,L+1); m_y = zeros(n,n-1,L+1); rho = zeros(n,n,L);

for l=0:L
    for i=1:n-1
        for j=1:n
            m_x(i,j,l+1) = u((j-1)*(n-1)+i+l*n*(n-1));
        end
    end
    for i=1:n
        for j=1:n-1
            m_y(i,j,l+1) = u((L+1)*(n-1)*n+(j-1)*n+i+l*n*(n-1));
        end
    end
end

for l = 1:L
    for i=1:n
        for j=1:n
            rho(i,j,l) = u(2*(L+1)*(n-1)*n+(j-1)*n+i+(l-1)*n*n);
        end
    end
end

end