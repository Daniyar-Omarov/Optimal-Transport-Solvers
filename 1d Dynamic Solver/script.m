clear; clc; Nt_list = [2; 4; 8]; N_list = [2; 4; 8; 16; 32; 64; 128];
cost = zeros(7,3); err = zeros(7,3); timeval = zeros(7,3);

cost_ex = 1/120; xa = 0; xb = 1; ya = 0; yb = 1; % Ex 1
% cost_ex = (25/16)*asin(4/5)-17/12; xa = 0; xb = 1; ya = 0; yb = 1; % Ex 2
% cost_ex = 1/30; xa = -1; xb = 1; ya = -1; yb = 1; % Ex 3
% cost_ex = 1.25; xa = -10; xb = 10; ya = -10; yb = 10; % Ex 4

for i=1:7
    for j=1:3
        tic; cost_val = lin_min(Nt_list(j), N_list(i), N_list(i), xa, xb, ya, yb);
        timeval(i,j) = toc; cost(i,j) = cost_val(2); err(i,j) = abs(cost(i,j) - cost_ex);
    end
end
