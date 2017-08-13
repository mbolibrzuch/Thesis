% This is a solver for differential equations. It uses the central difference
% method as well as (optionally) the forward difference method and compares
% the convergence of the two via loglog plots.
%
% 04-05: 
% This version tries to estimate the y' = sin(xy)+y, y(0)=1. So far I
% adjusted the code to:
%
% 1) Obtain the actual estimates, plot the best estimate, plot the 
% estimates for different step sizes. They look plausible, it's roughly the
% exponential function with a sinusoidal wave for lower values of x and y. 
% 
% 2) I save the vectors of estimates for different step sizes in order to
% be able to compare the values and plot convergence. The initial
% inspection of the matrix containing estimates indicates that things are
% most likely working.
% 
% Currently I'm trying to implement the actual comparison. I can't quite
% see which values to compare yet. By going in the way we discussed in the
% last meeting (1 to 1, 2 to 11, 3 to 21 etc.) the last value would be
% compared to the non-existing value (10 to 101, where the next vector in
% line contains 100 values). I will look more into how I'm actually
% indexing everything later to try and make it work.
%
% 04-06:
%
% I obtained an approximate for y(1) from another solver to help figure out
% which indexes don't make sense. It should be close to 3.51425. This leads
% me to believe that cutting the loop at i=n-2 is not correct, as values
% keep getting closer to the reference estimate for the i = n-1 and n. I
% will need to spend more time understanding what exactly my algorithm is
% doing.
% 
% 04-18:
%
% I figured out the indexing on paper. I'm attempting to implement the
% necessary changes and finally get this part to work.
%%
close all
clear all
n = 2;
max_m = 10;
m = 1;
% error = zeros(1,max_m);
save = zeros(n.^(max_m),max_m);
% error1 = zeros(1,max_m);

while m <= max_m
    n = 2.^m;  
    h = 1/n;
    H(m) = h;
    y = zeros(1,n+2);
%     x = zeros(1,n+2);
%     real = zeros(1,n+2);
    y(1) = 1;
%     y(2) = exp(h.^2/2);
% We can use the result of the worse algorithm with one step as a
% second initial value. 
% We can also use the derivative at the first point, in this case it's 1.
    y(2) = 1;
    save(1,m) = y(1);
    save(2,m) = y(2);
%     x(1) = 1;
%     x(2) = y(2);
%     real(1) = 1;
%     real(2) = y(2);
    for i = 1:n
%         x(i+2) = x(i+1) + h*(i+1)*h*x(i+1);
        y(i+2) = y(i) + 2*h*(sin(y(i+1)*h*(i+1))+y(i+1));
%         real(i+2) = exp((((i+2)*h).^2)/2);
    save(i+2,m) = y(i+2);
%     save1(i+2,m) = (i+1)*h;
    end   
    
% Errors
 
     
%     error1(m) = max(abs(x-real));
%     error(m) = max(abs(y-real));
    m = m+1;
end

% 04-18:
% There has got to be a better way to do this part, but for now I'm just
% going through the previous loops again recovering the differences I need.
% The loglog plot currently has a deformation from the straight line, but
% it does indicate convergence. I will further investigate the deformation
% later on and find out whether it is a problem in the code or a different
% result to interpret.

m = 1;
n = 2;

while m <= max_m - 1
n = 2.^m;
for i = 1:n
  error1(i) = abs(save(n+2,m)-save((2*n)+2,m+1));
end
error2(m) = max(error1);
m = m+1;
end

Hflip = fliplr(H(1:max_m-1));
errorflip = fliplr(error2);
plot(log2(Hflip), log2(errorflip))


% semilogy(error)
% hold on
% semilogy(error1)
% legend('Better', 'Worse');
% Check SUBPLOTS
% plot(y)
% plot(save)
% hold on
% plot(real)
% Check for y'(x) = sin(xy) + y
% y(0) = 1
% Plot against exp(x^2/2-x), exp(x)
% This below is for conventional form
% Hflip = fliplr(H);
% errorflip = fliplr(error);
% plot(log2(Hflip), log2(errorflip))
% Plot a straight line to see if its perfectly straight
