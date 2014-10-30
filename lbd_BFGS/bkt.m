function [alpha, count] = bkt(x,p,fcnHandle,alpha0, c, rho)
alpha = alpha0; count = 1;
[fcnval, gradval] = fcnHandle(x, 1);
temp = fcnHandle(x+alpha*p, 0);
while temp > fcnval + c * alpha * gradval' * p
	alpha = rho * alpha;
	temp = fcnHandle(x+alpha*p, 0);
	count = count + 1;
end
