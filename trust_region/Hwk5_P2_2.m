StopError = 1e-9; MaxIter = 1000;

TRCl_avg = zeros(2,1); TRCf_avg = zeros(2,1);
TRCl_std = zeros(2,1); TRCf_std = zeros(2,1);
LSCl_avg = zeros(2,1); LSCf_avg = zeros(2,1);
LSCl_std = zeros(2,1); LSCf_std = zeros(2,1);

for i = 1:2
	n = 200 * (i+3);
	TRCl = zeros(20,1); TRCf = zeros(20,1);
        LSCl = zeros(20,1); LSCf = zeros(20,1);
	for j = 1:20
		xstart = randn(n,1);
		[xTR, TRCl(j), TRCf(j)] = TrustRegion2(100, 0.1, StopError, MaxIter, xstart, @Rosenbrock_wrap, 0.5);
		[xLS, LSCl(j), LSCf(j)] = LineSearch2(1, 0.5, 0.5, StopError, MaxIter, xstart, @Rosenbrock_wrap, 0.5);
	end
	TRCl_avg(i) = mean(TRCl); TRCf_avg(i) = mean(TRCf);
	TRCl_std(i) = std(TRCl); TRCf_std(i) = std(TRCf);
	LSCl_avg(i) = mean(LSCl); LSCf_avg(i) = mean(LSCf);
        LSCl_std(i) = std(LSCl); LSCf_std(i) = std(LSCf);
end
