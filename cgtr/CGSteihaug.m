function [pd, k] = CGSteihaug(toler,g,B,Delta)
	z=zeros(length(g),1); r=g; d=-r; k=0;
	if norm(r)<toler
		pd = z; return
	end
	while 1>0
		k = k+1;
		tmp0 = B*d;
		tmp1 = dot(d, tmp0);
		if tmp1<=0
			c1=dot(d,d); c2=2*dot(d,z); c3=dot(z,z)-Delta^2;
			%sol = roots([c1 c2 c3]);
			%solreal = sol(imag(sol)==0);
			%o1=dot(g,d)+dot(tmp0,z); o2=tmp1/2;
			%model = o1*solreal + o2*solreal.^2;
			%[val, idx] = max(model)
			%tau = solreal(idx(1));
			tau = (-c2 + sqrt(c2*c2 - 4*c1*c3))/(2*c1);
			pd = z+tau*d; return
		end	
		tmp2 = dot(r,r);	
		alpha = tmp2/tmp1;
		z0 = z;
		z = z+alpha*d;
		if norm(z)>=Delta
			c1=dot(d,d); c2=2*dot(d,z0); c3=dot(z0,z0)-Delta^2;
			sol = roots([c1 c2 c3]);
 			tau = sol(imag(sol)==0 & sol>=0);
			pd = z0+tau(1)*d; return
		end
		r = r+alpha*tmp0; 
		if norm(r)<toler
			pd = z; return
		end
		beta = dot(r,r)/tmp2;
		d = -r+beta*d;
	end

