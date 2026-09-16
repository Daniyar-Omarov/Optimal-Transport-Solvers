clear; clc; Nt_list = [2; 4; 8]; N_list = [2; 4; 8; 16];
cost = zeros(7,3); err = zeros(7,3); timeval = zeros(7,3);

cost_ex = 0.013321449188603; xa = 0; xb = 1; % Ex 1
% cost_ex = 0.015511141991097; xa = -0.5; xb = 0.5; % Ex 2

for i=1:4
    for j=1:3
        tic; cost_val = lin_min_2d(Nt_list(j), N_list(i), xa, xb);
        timeval(i,j) = toc; cost(i,j) = cost_val(2); err(i,j) = abs(cost(i,j) - cost_ex);
    end
end