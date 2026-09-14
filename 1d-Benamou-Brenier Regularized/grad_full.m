function value = grad_full(n, L, beta, m, rho, rho0, dx)

% Defining Derivatives w.r.t m
value1 = zeros((n-1)*(L+1),1);
for l=0:L
    for i=1:n-1
        if l==0
            value1(l*(n-1)+i) = 4*m(i,l+1)/(rho0(i)+rho0(i+1));
        else
            value1(l*(n-1)+i) = 4*m(i,l+1)/(rho(i,l)+rho(i+1,l));
        end
    end
end

% Defining Derivativs w.r.t rho
value2 = zeros(L*n,1);
for l=1:L
    for i=1:n
        if i==1
            value2(i + (l-1)*n) = - 2*(m(i,l+1)^2)/(rho(i,l)+rho(i+1,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,l))-log(rho(i+1,l)))^2 + 2*(rho(i,l)+rho(i+1,l))*(log(rho(i,l))-log(rho(i+1,l)))/rho(i,l));
        elseif i==n
            value2(i + (l-1)*n) = - 2*(m(i-1,l+1)^2)/(rho(i-1,l)+rho(i,l))^2 + (beta^2/(2*dx^2))*((log(rho(i-1,l))-log(rho(i,l)))^2 - 2*(rho(i-1,l)+rho(i,l))*(log(rho(i-1,l))-log(rho(i,l)))/rho(i,l));
        else
            value2(i + (l-1)*n) = - 2*(m(i,l+1)^2)/(rho(i,l)+rho(i+1,l))^2 + (beta^2/(2*dx^2))*((log(rho(i,l))-log(rho(i+1,l)))^2 + 2*(rho(i,l)+rho(i+1,l))*(log(rho(i,l))-log(rho(i+1,l)))/rho(i,l)) - 2*(m(i-1,l+1)^2)/(rho(i-1,l)+rho(i,l))^2 + (beta^2/(2*dx^2))*((log(rho(i-1,l))-log(rho(i,l)))^2 - 2*(rho(i-1,l)+rho(i,l))*(log(rho(i-1,l))-log(rho(i,l)))/rho(i,l));
        end
    end
end

value = [value1; value2];

end