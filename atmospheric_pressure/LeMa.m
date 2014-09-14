function [x,k] = LeMa(StepBound, eta, xstart, StopError, method)
Delta = StepBound/2; x = xstart; E=Inf; k=0;
while E > StopError
	k = k+1;
	[f, g, B] = enso(x);
	switch method
	case 1
		[Q, S, Q] = svd(B);
		[lambda, pd] = EigenSolTR(Q, S, g, Delta, StopError);
		E = norm(pd);
	case 2
		[lambda, pd] = ExactTR(B, g, Delta, StopError);
		E = norm(pd);
	end
	x1 = x+pd;
	[f1, g1, B1] = enso(x1);
	rho = (f-f1) / (-dot(pd, g)-pd'*B*pd/2);
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
end 

