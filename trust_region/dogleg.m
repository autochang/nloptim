function [pd] = dogleg(B, g, Delta)
	pdB = -B \ g;
	pdU = -((g'*g) / (g'*B*g)) * g;
	pdiff = pdB - pdU;
	%if norm(pdB) < Delta
		pd = pdB;
	%elseif norm(pdU) > Delta
		pd = pdU / norm(pdU) * Delta;
	%else
	%	tmin1=roots([pdiff'*pdiff,2*pdiff'*pdU,pdU'*pdU - Delta*Delta]);
        %	tmin1 = tmin1(tmin1 >= 0 & tmin1 <= 1);
        %	assert(numel(tmin1) == 1);
        %	pd = pdU + (tmin1)*pdiff;	
	%end
	if norm(pdU) > Delta
		tau = Delta / norm(pdU);
	elseif norm(pdB) < Delta
		tau = 2;
	else
		r = roots([norm(pdB-pdU)^2 sum(pdB .* pdU *2, 1), norm(pdU)^2-Delta^2]);
		tau = r(r>=0) + 1;
	end
	if (tau >= 0 & tau <= 1)
		pd = tau * pdU;
	elseif (tau > 1 & tau <= 2)
		pd = pdU + (tau-1) * (pdB - pdU);
	else
		error('WRONG!')
	end

