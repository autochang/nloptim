function [G, c, l, u] = DataQPD(H, mu, lambda, C, Data)
	[n, m] = size(Data); m=m-1;
	G = H + mu*Data(:,1)*Data(:,1)';
	c = -ones(n,1) - lambda*Data(:,1);
	l = zeros(n,1);
	u = C*ones(n,1);
