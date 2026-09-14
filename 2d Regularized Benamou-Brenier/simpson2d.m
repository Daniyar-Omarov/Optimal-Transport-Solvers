function value = simpson2d(nd, a, b)

func = @(x1, x2) rho1_func(x1, x2)*log(rho1_func(x1,x2))-rho0_func(x1, x2)*log(rho0_func(x1,x2));
value = 0; x1 = linspace(a,b,nd+1); x2 = linspace(a,b,nd+1); h = 1/nd;

for i=1:nd
    g_low = 0; g_mid = 0; g_up = 0;
    for j=1:nd
        g_low = g_low + (h/6)*(func(x1(i), x2(j))+4*func(x1(i), x2(j)+(h/2))+func(x1(i), x2(j+1)));
        g_mid = g_mid + (h/6)*(func(x1(i)+(h/2), x2(j))+4*func(x1(i)+(h/2), x2(j)+(h/2))+func(x1(i)+(h/2), x2(j+1)));
        g_up = g_up + (h/6)*(func(x1(i+1), x2(j))+4*func(x1(i+1), x2(j)+(h/2))+func(x1(i+1), x2(j+1)));
    end
    value = value + (h/6)*(g_low + 4*g_mid + g_up);
end
end