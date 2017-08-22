% 08-22:
% 
% This version of the solver focuses on estimating the convergence of a
% problem with a known solution, specifically y' = xy, y(0)=1. It
% accomplishes the following:
% 1. Plots the convergence in the conventional, easy to read, logarithmic
% scale, for both the forward difference and central difference. Analytic
% calculations indicate that the first method should have the order of
% convergence one, and the second order of convergence two and the numerical results confirm it.
% 2. Checks the empirical order of convergence with polyfit MATLAB
% function.

close all
clear all
n = 2;
max_m = 13;
m = 1;
error = zeros(1,max_m);
error1 = zeros(1,max_m);

while m <= max_m
    n = 2.^m;  
    h = 1/n;
    H(m) = h;
    y = zeros(1,n+2);
    x = zeros(1,n+2);
    real = zeros(1,n+2);
    y(1) = 1;
    y(2) = 1;
    x(1) = 1;
    x(2) = y(2);
    real(1) = 1;
    real(2) = y(2);
    for i = 1:n-2
        x(i+2) = x(i+1) + h*(i+1)*h*x(i+1);
        y(i+2) = y(i) + 2*y(i+1)*h*h*(i+1);
        real(i+2) = exp((((i+2)*h).^2)/2);
    end
    error1(m) = max(abs(x-real));
    error(m) = max(abs(y-real));
    m = m+1;
end

Hflip = fliplr(H);
errorflip = fliplr(error);
errorflip1 = fliplr(error1);
plot(log2(Hflip), log2(errorflip))
hold on
plot(log2(Hflip), log2(errorflip1))
legend('Better' , 'Worse');
% % Check empirical convergence
Q = polyfit(log2(Hflip(1:12)), log2(errorflip(1:12)),1)
R = polyfit(log2(Hflip(1:12)), log2(errorflip1(1:12)),1)
