function [rho] = reductratio(x, pd, f, g, B, FHandle)
	[f1] = FHandle(x+pd,0);
	m = f;
	m1 = f + g'*pd + 0.5* pd'*B*pd;
	rho = (f-f1) / (m-m1);

