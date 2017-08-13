% 04-06:
%
% This version aims to use the solver to compute integrals. I use a simple
% example of sin(x). The solver allows me to solve the equation y' =
% sin(x), y(0) = k, with k being any coefficient, with relative ease. Since
% I just want to find the integral of sin(x) between 0 and 1, I will set
% k = -1. Then the solution to the equation is y = -cos(x) + c, with
% initial values we get: -1 = -cos(0) + c, c = 0. This way my algorithm
% will be computing the primitive of sin(x), namely -cos(x).
%
% 04-18:
% 
% This program needs a minor tweak so that it actually solves a particular
% problem. I want it to compute the integral of sin(x) from 0 to t, in this
% case I will set t to pi. This means I want to get the initial value back
% to 0 and in the end just substract the initial value.
%
% Something is not quite right with the error, however. Need to investigate
% further.
%%
close all
clear all
n = 2;
max_m = 17;
m = 1;
error = zeros(1,max_m);
error1 = zeros(1,max_m);

while m <= max_m
    n = 2.^m;  
    h = pi/n;
    H(m) = h;
    y = zeros(1,n+2);
    x = zeros(1,n+2);
    real = zeros(1,n+2);
    y(1) = 0;
%     y(2) = exp(h.^2/2);
% We can use the result of the worse algorithm with one step as a
% second initial value. We can also use the derivative at the initial value
% to estimate the second value.
    y(2) = 0;
    x(1) = 0;
    x(2) = y(2);
    real(1) =  0;
    real(2) = y(2);
    for i = 1:n-2
        x(i+2) = x(i+1) + sin(h*(i+1))*h;
        y(i+2) = y(i) + 2*h*sin((i+1)*h);
        real(i+2) = -cos((i+2)*h)+1;
    end
    
    if m >= 2
    error1(m-1) = abs(x(n)-real(n));
    error(m-1) = abs(y(n)-real(n));
    end
    m = m+1;
    
end
% 
% semilogy(error)
% hold on
% semilogy(error1)
% legend('Better', 'Worse');
% Check SUBPLOTS
% plot(y(1:n))
% hold on
% plot(real(1:n))
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
% integral = y(n)-y(1)
