function [xsol, NumSV] = BCAugLag(xstart, lambda0, etaTol, omegaTol, C, Data)
	mu = 10; omega = 1/mu; eta = 1/(mu^(0.1)); NumIter = 0;
	x = xstart; lambda = lambda0; n=length(x);
	y = Data(:,1); H = MatQPD(Data);
	while 1>0
		NumIter = NumIter + 1;
		% generate the model coef for the updated mu and lambda
		[G, c, l, u] = DataQPD(H, mu, lambda, C, Data);
		% find approx sol of the subproblem 17.50 satisfying ||...|| <= omega_k
		[x,fct,ier,nsub]=minq(0,c,G,l,u,[],ones(n,1));
		gd = G*x+c;
		% update lambda and mu
		if abs(dot(x,y)) <= eta
			stop1 = (abs(dot(x,y)) <= etaTol);
			Proj = ProjBox(x-gd, l, u);
			stop2 = (norm(x-Proj) <= omegaTol);
			if stop1 & stop2
				xsol = x;
				NumSV = sum(abs(xsol)>1e-9); 
				return
			end
			lambda = lambda - mu*dot(x,y);
			mu = mu;
			eta = eta / (mu)^(0.9);
			omega = omega / mu;
		else
			lambda = lambda;
			mu = 100*mu;
			eta = 1 / (mu)^(0.1);
			omega = 1 / mu;
		end
	end
