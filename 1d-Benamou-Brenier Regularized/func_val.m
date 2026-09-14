function value = func_val(n, L, beta, m, rho, rho0, dx)

value = 0; 
for l=0:L
    for i=1:n-1
        if l==0
            value = value + 2*(m(i,l+1)^2)/(rho0(i)+rho0(i+1)) + ((log(rho0(i))-log(rho0(i+1)))^2)*(rho0(i)+rho0(i+1))*(beta^2)/(2*dx^2);
        else
            value = value + 2*(m(i,l+1)^2)/(rho(i,l)+rho(i+1,l)) + ((log(rho(i,l))-log(rho(i+1,l)))^2)*(rho(i,l)+rho(i+1,l))*(beta^2)/(2*dx^2);
        end
    end
end

end