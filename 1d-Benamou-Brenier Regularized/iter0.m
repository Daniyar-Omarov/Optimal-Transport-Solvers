function [m, rho] = iter0(rho0, rho1, n, L, dx)

% n - grid size; L+1 - number of interval in time;
% dx - mesh size in x and y; dt - mesh size in time;
dt = 1/(L+1); t_l = dt*linspace(0, L+1, L+2);

rho = zeros(n,L);
for l=1:L
    for i=1:n
        rho(i,l) = (1-t_l(l+1))*rho0(i) + t_l(l+1)*rho1(i);
    end
end

m = zeros(n-1,L+1);
for l=0:L
    if l==0
        m(1,l+1) = - dx*(rho(1,l+1)-rho0(1))/(dt);
    elseif l==L
        m(1,l+1) = - dx*(rho1(1)-rho(1,l))/(dt);
    else
        m(1,l+1) = - dx*(rho(1,l+1)-rho(1,l))/(dt);
    end
    for i=2:n-1
        if l==0
            m(i,l+1) = m(i-1,l+1) - dx*(rho(i,l+1)-rho0(i))/(dt);
        elseif l==L
            m(i,l+1) = m(i-1,l+1) - dx*(rho1(i)-rho(i,l))/(dt);
        else
            m(i,l+1) = m(i-1,l+1) - dx*(rho(i,l+1)-rho(i,l))/(dt);
        end
    end
end

end