function toss = randi_expn(tau,n)

toss = zeros(n,1);

x_max = 100;
xx = -x_max:1:x_max;
yy = 1000*exp(-(0:x_max)/tau);
yy = [yy(end:-1:2) yy];
rr = randi([0 floor(sum(yy))],n);

for ii = 1:n
    toss(ii) = xx(find(cumsum(yy)>=rr(ii),1,'first'));
end

%for testing
%nn = histcounts(toss,xx-0.5);
%plot(xx(1:end-1),nn,'k')