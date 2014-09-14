function xx = plp(t, x, gd, tbar)
	n = length(x); xx=zeros(n,1);
	for i=1:n
		if t <= tbar(i)
			xx(i,1) = x(i)-t*gd(i);
		else
			xx(i,1) = x(i)-tbar(i)*gd(i);
		end
	end
