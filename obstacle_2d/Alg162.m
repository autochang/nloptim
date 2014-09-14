function [xplus, count] = Alg162(xstart, G, c, P, l, u, toler)
	x = xstart; count = 0;
	r = G*x+c; g = P*r; d = -g;
	rgdot = Inf;
	while 1>0
		count = count + 1;
		Gd = G*d;
		alpha = dot(r,g)/dot(d,Gd);
		x = x+alpha*d;
		stopper = prod(x>=l)*prod(x<=u);
		rplus = r+alpha*Gd;
		gplus = P*rplus;
		beta = dot(rplus,gplus)/dot(r,g);
		d = -gplus+beta*d;
		g = gplus; r = rplus; rgdot = dot(r,g);
		if stopper == 1 | rgdot < toler
                        xplus = x; return
                end 
	end
