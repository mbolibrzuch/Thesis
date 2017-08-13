% 10-05
%
% This version attempts to show the problems with having a non-smooth
% function within the problem. We use (x^2)^(1/3) as f(x) and solve the DE:
% y' = f(x), y(0) = 1.
%
% I tried a few versions, with real value as well as with convergence
% between approximation and I could not obtain anything other than
% convergence expected for a smooth enough function.
%
% I suspect it could be because of the fact that I'm only looking for a
% solution on the (0,1) interval and the cusp of this function is at 0. I
% tried modifying the function to have the cusp at 1/sqrt(2) instead but the
% error is still similar, so perhaps that is not the problem.
%
% I'm considering the possibility that the central difference method does
% not require fully smooth functions in practical implementation. 
% Then I would need to:
% 1) Understand why the non-smoothness does not affect this method
% 2) Find a method that does require smooth functions and understand what
% happens with the error there
% 
% However, that seems rather odd, considering that the Taylor polynomial
% requires a function to be differentiable to exist. A function with a cusp
% is not differentiable at all at the cusp so it should exhibit some
% issues.
%
% 05-11:
% 
% After looking at this on the following day I thought that perhaps I want
% to just take a derivative of the chosen f(x) using central difference
% method. Does not work so far, I will proceed to do more calculations on
% paper to get a better idea of what I want to happen. 

close all
clear all
n = 2;
max_m = 15;
m = 1;
% error = zeros(1,max_m);
save = zeros(n.^(max_m),max_m);

while m <= max_m
    n = 2.^m;
    h = 1/n;
    H(m) = h;
    y = zeros(1,n+2);
    %     real = zeros(1,n+2);
%     y(1) = 1;
%     y(2) = 1;
    %     x(1) = 1;
    %     x(2) = y(2);
    save(1,m) = y(1);
    save(2,m) = y(2);
    %     real(1) = 1;
    %     real(2) = y(2);
    
    for i = 1:n
        %%
        % These are the different forms of inputs I have tried so far.
        % 1. Cusp at 1/sqrt(2):
%                 y(i+2) = y(i) + 2*h*((h*(i+1))-(1/(sqrt(2))).^2).^(1/3);
        % 2. |x-1/sqrt(2)|, like in the KTH document to try and mimic it:
                y(i+2) = y(i) + 2*h*abs((h*(i+1))-(1/sqrt(2)));
        % 3. |x-1/sqrt(2)|, just computing a derivative. I believe this is
        % wrong for now.
%                 y(i+1) = (abs((i+2)*h)-(1/sqrt(2))-((i*h)-1/sqrt(2)))/(2*h);
        % I tried both real values and error between approximation.
        %   real(i+2) = ( (3/5)*(((h*(i+1)).^2).^(1/3))*h*(i+1) ) + 1;
        save(i+2,m) = y(i+2);
        
    end
    
    % error1(m) = max(abs(x-real));
    %     error(m) = max(abs(y-real));
    m = m+1;
end

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
% plot(y(1:n-2))
% hold on
% plot(real)
% Check for y'(x) = sin(xy) + y
% y(0) = 1
% Plot against exp(x^2/2-x), exp(x)
% This below is for conventional form
% Hflip = fliplr(H);
% errorflip = fliplr(error);
% errorflip1 = fliplr(error1);
% plot(log2(Hflip), log2(errorflip))
% hold on
% plot(log2(Hflip), log2(errorflip1))
% legend('Better' , 'Worse');
% % Plot a straight line to see if its perfectly straight
% Q = polyfit(log2(Hflip(1:12)), log2(errorflip(1:12)),1)
% R = polyfit(log2(Hflip(1:12)), log2(errorflip1(1:12)),1)