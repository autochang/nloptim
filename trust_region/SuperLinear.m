function [ErrorQuotient] = SuperLinear(x, xout)
	TotalNum = size(x,2);
	ErrorQuotient = zeros(TotalNum-1, 1);
	Error = zeros(TotalNum, 1);
	Diff = x - repmat(xout, 1, TotalNum);
	Err1 = ( sqrt( sum( Diff(:,1:TotalNum-1).^2 ) ) )';
	Err2 = ( sqrt( sum( Diff(:,2:TotalNum).^2 ) ) )';
	%LogErrorQuotient = log10(Err2) - log10(Err1);
	I1 = (1:(TotalNum-1))'; I2 = (2:TotalNum)';
	K1 = I1 .^ I1; K2 = I2 .^ I2;
	ErrorQuotient = (Err2 .* K2) ./ (Err1 .* K1);
	%ErrorQuotient = Err2 ./ Err1;
	
