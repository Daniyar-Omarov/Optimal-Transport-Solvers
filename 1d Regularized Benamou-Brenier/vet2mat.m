function [m, rho] = vet2mat(n, L, u)

m = zeros(n-1,L+1); rho = zeros(n,L);
for l = 1:L+1
    for i=1:n-1
        m(i,l) = u(i+(l-1)*(n-1));
    end
end
for l=1:L
    for i=1:n
        rho(i,l) = u((L+1)*(n-1)+i+(l-1)*n);
    end
end

end