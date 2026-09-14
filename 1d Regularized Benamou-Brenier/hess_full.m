function value = hess_full(n, L, beta, m, rho, rho0)

% Defining the Second Derivative w.r.t. m
value1 = zeros((n-1)*(L+1),1);
for l=0:L
    for i=1:n-1
        if l==0
            value1(i) = 4/(rho0(i)+rho0(i+1));
        else
            value1(l*(n-1)+i) = 4/(rho(i,l)+rho(i+1,l));
        end
    end
end

% Defining Mixed Derivative w.r.t. m and rho
value2 = zeros((n-1)*(L+1),L*n);
for l = 1:L
    for i=1:n-1
            value2((n-1)+i+(l-1)*(n-1), i+(l-1)*n) = -4*m(i,l+1)/(rho(i,l)+rho(i+1,l))^2;
            value2((n-1)+i+(l-1)*(n-1), i+1+(l-1)*n) = -4*m(i,l+1)/(rho(i,l)+rho(i+1,l))^2;
    end
end

value3 = zeros(L*n);
for l = 1:L
    value3(1+(l-1)*n:l*n,1+(l-1)*n:l*n) = J_rr(rho(:,l), m(:,l+1), n) + (beta^2)*I_rr(rho(:,l), n);
end

value = zeros((n-1)*(L+1) + n*L);
value(1:(n-1)*(L+1),(n-1)*(L+1)+1:end) = value2;
value = value+value';
for i=1:(n-1)*(L+1)
    value(i,i) = value1(i);
end
value((n-1)*(L+1)+1:end,(n-1)*(L+1)+1:end) = value3;

end