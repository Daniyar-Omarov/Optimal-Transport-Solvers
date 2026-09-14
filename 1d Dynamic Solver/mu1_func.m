function value = mu1_func(x)

% Example 1
value = (4 * x) .* (x < 0.5) + (4 * (1 - x)) .* (x >= 0.5);

% Example 2
value = (2*x+1)/2;

% Example 3
value = (1/2)*ones(size(x));

% Example 4
value = normpdf(x, 0, 1);

end