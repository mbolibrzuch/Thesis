
m = 1;
max_m1 = 11;
n = 2;
while m <= max_m1 - 1
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


