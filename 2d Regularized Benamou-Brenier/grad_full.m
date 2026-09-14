function value = grad_full(n, L, dx, beta, m_x, m_y, rho, rho0)

% Defining Derivativs w.r.t m^x and m^y
value1 = zeros(2*n*(n-1)*(L+1),1);
for l=0:L
    if l==0
        for j=1:n
            for i=1:n-1
                value1((j-1)*(n-1)+i) = 4*m_x(i,j,l+1)/(rho0(i,j)+rho0(i+1,j));
            end
        end
        for i=1:n
            for j=1:n-1
                value1((L+1)*n*(n-1)+(j-1)*n+i) = 4*m_y(i,j,l+1)/(rho0(i,j)+rho0(i,j+1));
            end
        end
    else
        for i=1:n-1
            for j=1:n
                value1(l*(n-1)*n+(j-1)*(n-1)+i) = 4*m_x(i,j,l+1)/(rho(i,j,l)+rho(i+1,j,l));
            end
        end
        for i=1:n
            for j=1:n-1
                value1((L+1)*n*(n-1)+l*(n-1)*n+(j-1)*n+i) = 4*m_y(i,j,l+1)/(rho(i,j,l)+rho(i,j+1,l));
            end
        end
    end
    
end

% Defining Derivativs w.r.t rho
value2 = zeros(L*n^2,1);
for l=1:L
    for j=1:n
        for i=1:n
            if i==1
                value2(i + (j-1)*n + (l-1)*n*n) = - 2*(m_x(i,j,l+1)^2)/(rho(i,j,l)+rho(i+1,j,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,j,l))-log(rho(i+1,j,l)))^2 + 2*(rho(i,j,l)+rho(i+1,j,l))*(log(rho(i,j,l))-log(rho(i+1,j,l)))/rho(i,j,l));
            elseif i==n
                value2(i + (j-1)*n + (l-1)*n*n) = - 2*(m_x(i-1,j,l+1)^2)/(rho(i-1,j,l)+rho(i,j,l))^2 + (beta^2/(2*dx^2))*((log(rho(i-1,j,l))-log(rho(i,j,l)))^2 - 2*(rho(i-1,j,l)+rho(i,j,l))*(log(rho(i-1,j,l))-log(rho(i,j,l)))/rho(i,j,l));
            else
                value2(i + (j-1)*n + (l-1)*n*n) = - 2*(m_x(i,j,l+1)^2)/(rho(i,j,l)+rho(i+1,j,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,j,l))-log(rho(i+1,j,l)))^2 + 2*(rho(i,j,l)+rho(i+1,j,l))*(log(rho(i,j,l))-log(rho(i+1,j,l)))/rho(i,j,l)) - 2*(m_x(i-1,j,l+1)^2)/(rho(i-1,j,l)+rho(i,j,l))^2 + (beta^2/(2*dx^2))*((log(rho(i-1,j,l))-log(rho(i,j,l)))^2 - 2*(rho(i-1,j,l)+rho(i,j,l))*(log(rho(i-1,j,l))-log(rho(i,j,l)))/rho(i,j,l));
            end
        end
    end
end

value3 = zeros(L*n^2,1);
for l=1:L
    for j=1:n
        for i=1:n
            if j==1
                value3(i + (j-1)*n + (l-1)*n*n) = - 2*(m_y(i,j,l+1)^2)/(rho(i,j,l)+rho(i,j+1,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,j,l))-log(rho(i,j+1,l)))^2 + 2*(rho(i,j,l)+rho(i,j+1,l))*(log(rho(i,j,l))-log(rho(i,j+1,l)))/rho(i,j,l));
            elseif j==n
                value3(i + (j-1)*n + (l-1)*n*n) = - 2*(m_y(i,j-1,l+1)^2)/(rho(i,j-1,l)+rho(i,j,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,j-1,l))-log(rho(i,j,l)))^2 - 2*(rho(i,j-1,l)+rho(i,j,l))*(log(rho(i,j-1,l))-log(rho(i,j,l)))/rho(i,j,l));
            else
                value3(i + (j-1)*n + (l-1)*n*n) = - 2*(m_y(i,j,l+1)^2)/(rho(i,j,l)+rho(i,j+1,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,j,l))-log(rho(i,j+1,l)))^2 + 2*(rho(i,j,l)+rho(i,j+1,l))*(log(rho(i,j,l))-log(rho(i,j+1,l)))/rho(i,j,l)) - 2*(m_y(i,j-1,l+1)^2)/(rho(i,j-1,l)+rho(i,j,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,j-1,l))-log(rho(i,j,l)))^2 - 2*(rho(i,j-1,l)+rho(i,j,l))*(log(rho(i,j-1,l))-log(rho(i,j,l)))/rho(i,j,l));
            end
        end
    end
end
value4 = value2 + value3; value = [value1; value4];

end