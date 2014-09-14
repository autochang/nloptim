function [lambda, pd] = ExactTR(B, g, Delta, StopError)
	n = length(g); pd = B \ (-g); lambda = 0; E=Inf;
	if norm(pd) <= Delta
		return
	else
		while E > StopError
			lambda0 = lambda;
			R = chol(B+lambda*eye(n));
			p = -R\(R'\g);
			q = R'\p;
			lambda = lambda + (norm(p)/norm(q))^2 * ((norm(p)-Delta)/Delta);
			E = abs(lambda-lambda0);	
		end
		pd = p;
	end

