%% Plotting convergence of error in trapezoidal method.
% Goal: show that convergence is NOT of order 2 (expected order) when 
% computing for a function that has a cusp and making the grid points
% on-and off the cusp (big error-small error oscillation). Using built-in
% trapz function.
%%
m = 1;
max_m = 20;
error = zeros(1,m);
while m <= max_m
    n = 2.^m;
    h = 1/n;
    H(m) = h;
    X = 0:h:1;
    Y = ((X- 1/sqrt(2)).^2).^(1/3);
    value(m) = trapz(X,Y);
    if m >= 2
    error(m) = abs(value(m)-value(m-1));
    end
    m = m+1;
end

%% Plotting the convergence of error in conventional log-log form
%%
Hflip = fliplr(H(1:max_m-1));
errorflip = fliplr(error(2:max_m));
plot(log2(Hflip), log2(errorflip))