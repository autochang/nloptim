function y = Rosenbrock(x)
	n = length(x);
	xodd = x(1:2:n);
	xeven = x(2:2:n);
	y = sum( (ones(n/2, 1)-xodd).^2 + 10 * (xeven - xodd.^2).^2 );
