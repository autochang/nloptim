function [x, Cl, Cf] = TrustRegion(StepBound, eta, StopError, MaxIter, xstart, FHandle, delta)
k = 0; E = Inf; Delta = StepBound/2;
x = xstart; B = zeros(length(x), length(x));
Cl = 0; Cf = 0;
while (k < MaxIter) & (E > StopError)
	x0 = x;
	[f, g, H] = FHandle(x, 2); % one function eval
	[L, D, p] = ldl(H, 'vector');
	[Q, S, Q] = svd(full(D));
	F = mdf_sparse(Q, diag(S), delta);
	B(p,p) = L * (D+F) * L'; 
	pd = dogleg(B, g, Delta); % one linear solve
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
	E = norm(x-x0);
	k = k + 1;
	Cl = Cl + 1; Cf = Cf + 2;
end
