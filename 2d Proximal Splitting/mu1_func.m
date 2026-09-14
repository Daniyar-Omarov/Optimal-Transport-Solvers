function value = mu1_func(x1,x2)

% Example 1
% value = 1.631545990171433*exp(-2*(x1-0.25).^2-2*(x2-0.75).^2);

% Example 2
a0 = -1; a1 = 1; a2 = -a0*exp(-1/8);
b0 = -1; b1 = 1; b2 = -b0*exp(-1/8);
value = (a1 + a2*x1.*exp(0.5*x1.^2) + 0.01*pi*sin(pi*x1).*sin(pi*x2)).*(b1 + b2*x2.*exp(0.5*x2.^2) + 0.01*pi*sin(pi*x1).*sin(pi*x2)) - (0.01*pi*cos(pi*x1).*cos(pi*x2)).^2;

end