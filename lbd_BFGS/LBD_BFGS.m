function [x, k] = LBD_BFGS(xstart, m, TolError, fcnHandle)
	x=xstart; k=0; E=Inf; alpha0=1; c=0.5; rho=0.5;
	n=length(xstart);
	S=zeros(n,m); Y=zeros(n,m);
	Theta = zeros(m,1); Yhat=zeros(n,m); Bs=zeros(n,m);
	while E > TolError
		[f0, g0] = fcnHandle(x, 1); 
		if k < m
			if k == 0
				H0 = eye(n);
			else
				H0 = dot(S(:,k),Y(:,k)) / dot(Y(:,k),Y(:,k)) * eye(n);
			end
			B0 = diag(1./ diag(H0));
			if k == 0
				pd = -g0;
			else 
				Stmp = S(:, 1:k); Yhattmp = Yhat(:, 1:k);
				pd = LimBFGS(g0, Stmp, Yhattmp, H0);
			end
			[alpha, count] = bkt(x,pd,fcnHandle,alpha0,c,rho);
			x = x + alpha * pd;
			[f1, g1] = fcnHandle(x, 1);
			s = alpha*pd; y=g1-g0;
			E = norm(g1);
			S(:, k+1) = s;
			Y(:, k+1) = y;
			if k == 0
				if dot(s,y) >= 0.2*dot(s,s)
					Theta(k+1)=1;
				else
					Theta(k+1)=0.8*dot(s,s)/ (dot(s,s)-dot(s,y));
				end
				Bs(:,k+1) = s;
				Yhat(:,k+1)=Theta(k+1)*Y(:,k+1)+(1-Theta(k+1))*s;
			else
				Stmp = S(:, 1:k); Yhattmp = Yhat(:, 1:k); Bstmp = Bs(:, 1:k); Thetatmp = Theta(1:k);
				[Bs(:,k+1), Theta(k+1), Yhat(:,k+1)] = DamBFGS(Stmp, Bstmp, Thetatmp, Yhattmp, s, y, B0);
			end
			clearvars Stmp Ytmp Yhattmp Thetatmp;
		else
			H0 = dot(S(:,end),Y(:,end)) / dot(Y(:,end),Y(:,end)) * eye(n); % 7.20
			B0 = diag(1./ diag(H0));
			pd = LimBFGS(g0, S, Yhat, H0);
			[alpha, count] = bkt(x,pd,fcnHandle,alpha0,c,rho); % backtracking
			x = x + alpha * pd;
			[f1, g1] = fcnHandle(x, 1);
			E = norm(g1);
			s = alpha * pd; y = g1 - g0;
			[Dampedbs, Dampedtheta, DampedY] = DamBFGS(S, Bs, Theta, Yhat, s, y, B0);
			% update the new s and y
			S(:, 1:m-1) = S(:, 2:end);
			Y(:, 1:m-1) = S(:, 2:end);
			Theta(1:m-1) = Theta(2:end);
			Yhat(:, 1:m-1) = Yhat(:, 2:end);
			Bs(:, 1:m-1) = Bs(:, 2:end);
			S(:, end) = s;
			Y(:, end) = y;
			Theta(end) = Dampedtheta;
			Yhat(:, end) = DampedY;
			Bs(:, end) = Dampedbs;
		end 
		k = k+1;
	end

