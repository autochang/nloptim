function [xsol, k1, k2, k3] = GPCG(xstart, G, c, l, u, method, precond, toler)
	x = xstart; k1=0; k2=0; k3=0; breakpoint=Inf; gradient=Inf;
	n = length(x);
	while 1>0
		k3 = k3+1; % record the number of GP iterations
		gd = G*x+c;
		if abs(gradient)<=toler | breakpoint==0 % stopping rule given in the lecture
			xsol = x;
			return
		end
		[xcauchy, countcauchy, breakpoint, gradient]= FindCauchy(gd,x,l,u,G,c);% find the Cauchy point xc
		k1 = k1 + countcauchy; % record the number of Cauchy points iterations
		active = find(xcauchy==l | xcauchy==u);
		inactive = find(xcauchy~=l & xcauchy~=u);
		if isempty(active)+isempty(inactive) == 0
			Idm = eye(n,n);
			A = Idm(:,active)'; Z = Idm(:,inactive); b=xcauchy(active,1);
			switch precond
				case 0 
					nn = length(inactive);
					W = eye(nn);
				case 1 
					H = eye(n); % actually the same as case 0 due to Z
					W = Z'*H*Z;
				case 2
					H = diag(diag(abs(G)));
					W = Z'*H*Z;
				case 3
					H = G;
					W = Z'*G*Z;
			end
			P = Z*(W\Z'); % projection matrix in Algorithm 16.2
			% find approx sol x+ of (16.74) such that q(x+)<=q(xc) and x+ is feasible
			switch method
				case 0 % without CG part
					xplus = xcauchy; countCG = 0;
				case 1 % CG: Algorithm 16.1
					[xplus, countCG] = Alg161(xcauchy, A, b, Z, G, c, l, u, W, toler);
				case 2 % CG: Algorithm 16.2
					[xplus, countCG] = Alg162(xcauchy, G, c, P, l, u, toler);
				otherwise
					disp('Specify your method to solve the subproblem.')
			end
			k2 = k2 + countCG; % record the number of CG iterations
			x = xplus;
		else
			x = xcauchy;
		end
	end
