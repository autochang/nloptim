function [x, k, kk, e] = TR_HWK8(StepBound, eta, StopError, xstart, FHandle)
k = 0; E = Inf; Delta = StepBound/2; k=0; kk=0;
x = xstart; n=length(x); B = zeros(n,n);
while E > StopError
	[f, g, B] = FHandle(x, 2); % one function eval
	toler = min(0.5, sqrt(norm(g)))*norm(g);
	[pd, CGk] = CGSteihaug(toler,g,B,Delta);
	E = norm(g);
        k = k + 1;
        kk = kk + CGk;
	if pd == zeros(n,1)
		x = x; return
	end
	rho = reductratio(x, pd, f, g, B, FHandle); % one function eval
	if rho < 0.25
		Delta = Delta / 4;
	elseif (rho > 0.75) & (norm(pd) == Delta)
		Delta = min(2*Delta, StepBound);
	else
		Delta = Delta;
	end
	if rho > eta
		x = x + pd;
	else
		x = x;
	end
	e(k) = norm(x-ones(n,1));
end
