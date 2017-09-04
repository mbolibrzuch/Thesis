n = 1;
max_n = 5;
figure

while n <= max_n
    if mod(n,2) == 0
        X = [0-1/n 0 0+1/n];
        Y = [0 1 0];
        subplot(max_n,1,n)
        plot(X,Y);
        title(['T_' num2str(n)'])
        hold on
    else
      X = [1-1/n 1 1+1/n];
        Y = [0 1 0];
        subplot(max_n,1,n)
        plot(X,Y);
        title(['T_' num2str(n)'])
        hold on
    end
    n = n+1;
end

