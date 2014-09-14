function [x, k, e] = LeMaQN(StepBound, eta, xstart, StopError, method, xsol)
Delta = StepBound/2; x = xstart; E=Inf; k=0; n=length(x); QNS = 0.1*eye(n); 
while E > StopError
	x0 = x;
	k = k+1;
	[f, g, Resi, Jaco] = ensoQN(x);
	B = Jaco'*Jaco + QNS;
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
	[f1, g1, Resi1, Jaco1] = ensoQN(x1);
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
	e(k) = norm(x-xsol);
	if x == x0
		QNS = QNS;
	else
		s = x - x0;
		y = g1 - g;
		ypound = g1 - Jaco'*Resi1;
		tmp4 = abs(dot(s,ypound)/(s'*QNS*s));
                tau = min(1, tmp4);
                QNS = tau*QNS;
		tmp1 = ypound - QNS*s;
		tmp2 = dot(y,s);
		tmp3 = y*y';
		comp1 = (tmp1*y'+y*tmp1') / tmp2;
		comp2 = dot(tmp1, s)/(tmp2^2) * tmp3;
		QNS = QNS + comp1 - comp2;
	end
end 

