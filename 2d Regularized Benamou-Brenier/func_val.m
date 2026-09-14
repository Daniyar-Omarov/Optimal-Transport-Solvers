function value = func_val(n, L, dx, beta, m_x, m_y, rho, rho0)

value = 0;

for l=0:L
    
    if l==0
        for j=1:n
            for i=1:n-1
                value = value + 2*(m_x(i,j,l+1)^2)/(rho0(i,j)+rho0(i+1,j)) + ((log(rho0(i,j))-log(rho0(i+1,j)))^2)*(rho0(i,j)+rho0(i+1,j))*(beta^2)/(2*dx^2);
            end
        end
        
        for i=1:n
            for j=1:n-1
                value = value + 2*(m_y(i,j,l+1)^2)/(rho0(i,j)+rho0(i,j+1)) + ((log(rho0(i,j))-log(rho0(i,j+1)))^2)*(rho0(i,j)+rho0(i,j+1))*(beta^2)/(2*dx^2);
            end
        end
    else
        for j=1:n
            for i=1:n-1
                value = value + 2*(m_x(i,j,l+1)^2)/(rho(i,j,l)+rho(i+1,j,l)) + ((log(rho(i,j,l))-log(rho(i+1,j,l)))^2)*(rho(i,j,l)+rho(i+1,j,l))*(beta^2)/(2*dx^2);
            end
        end
        
        for i=1:n
            for j=1:n-1
                value = value + 2*(m_y(i,j,l+1)^2)/(rho(i,j,l)+rho(i,j+1,l)) + ((log(rho(i,j,l))-log(rho(i,j+1,l)))^2)*(rho(i,j,l)+rho(i,j+1,l))*(beta^2)/(2*dx^2);
            end
        end
        
    end
    
end

end