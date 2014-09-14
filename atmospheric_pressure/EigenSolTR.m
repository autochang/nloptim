function [lambda, pd] = EigenSolTR(Q, S, g, Delta, StopError)
	SS = diag(S); n=length(SS); E = Inf;
	lambda1 = min(SS);
	lambda0 = -lambda1 + rand(1);
	c = 0.5;
	lambda = lambda0; 
	linear0 = Q'*g ./ SS;
	if prod(SS) ~= 0 & norm(linear0) <= Delta
		lambda = 0;
	  	pd = -Q*linear0; 
	else
		while E > StopError
			lambdapre = lambda;
			linear1 = Q'*g ./ (SS+lambda*ones(n,1));
			phi = 1/Delta - 1/norm(linear1);
			linear2 = (Q'*g).^2 ./ (SS+lambda*ones(n,1)).^3; 
			temp = sum(linear2);
			phid = -temp / norm(linear1)^3;
			lambdatmp = lambda - phi / phid;
				if lambdatmp > -lambda1
					lambda = lambdatmp;
				else
					lambda = -c*lambda1+(1-c)*lambda;
				end
			E = norm(lambda-lambdapre);
		end
		linear = Q'*g ./ (SS+lambda*ones(n,1));
		pd = -Q*linear; 
	end
