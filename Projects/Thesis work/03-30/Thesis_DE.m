% 03-30:
% 
% This version of the solver focuses on estimating the convergence of a
% problem with a known solution, specifically y' = xy, y(0)=1. It
% accomplishes the following:
% 1. Plots the estimated solution of the equation.
% 2. Plots the convergence in the conventional, easy to read, logarithmic
% scale, for both the (forward difference NOTE: I'm still unsure about
% what to call the methods exactly) and (central difference). Analytic
% calculations indicate that the first method should have the order of
% convergence one, and the second order of convergence two.
% 3. Checks the empirical order of convergence with polyfit MATLAB
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
%     y(2) = exp(h.^2/2);
% We can use the result of the worse algorithm with one step as a
% second initial value. We can also use the derivative at the initial value
% to estimate the second value.
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
% 
% semilogy(error)
% hold on
% semilogy(error1)
% legend('Better', 'Worse');
% Check SUBPLOTS
% plot(y(1:n-2))
% hold on
% plot(real)
% Check for y'(x) = sin(xy) + y
% y(0) = 1
% Plot against exp(x^2/2-x), exp(x)
% This below is for conventional form
Hflip = fliplr(H);
errorflip = fliplr(error);
errorflip1 = fliplr(error1);
plot(log2(Hflip), log2(errorflip))
hold on
plot(log2(Hflip), log2(errorflip1))
legend('Better' , 'Worse');
% % Plot a straight line to see if its perfectly straight
Q = polyfit(log2(Hflip(1:12)), log2(errorflip(1:12)),1)
R = polyfit(log2(Hflip(1:12)), log2(errorflip1(1:12)),1)
