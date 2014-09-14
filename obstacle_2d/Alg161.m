function [xplus,count]=Alg161(xcauchy, A, b, Z, G, c, l, u, W, toler)
	count = 0;
	Yxy = A'*b;
	xz = Z'*(xcauchy-Yxy);
	xstar = Yxy + Z*xz;
	cz = Z'*G*Yxy + Z'*c;
	rz = Z'*G*Z*xz + cz;
	gz = W \ rz;
	dz = -gz;
	rquad = Inf;
	while 1>0
		count = count + 1;
		Zdz = Z*dz;
		GZdz = G*Zdz;
		alpha = dot(rz, gz) / dot(Zdz, GZdz);
		xz = xz + alpha*dz;
		xstar = Yxy + Z*xz;
		stopper = prod(xstar>=l)*prod(xstar<=u);
		rzplus = rz + alpha*Z'*GZdz;
		gzplus = W \ rzplus;
		beta = dot(rzplus,gzplus)/dot(rz,gz);
		dz = -gzplus + beta*dz;
		gz = gzplus; rz = rzplus; rquad = dot(rz, gz);
		if stopper == 1 | rquad < toler
                        xplus = xstar; return
                end 
	end
