StepBound = 10.^[-2;-1;0;1;2];
eta = 0.2;
xout = newtonLikeMethod(@cute_wrap, ones(100,1), 1, 1e-9, 1000);
for i = 1:5
	[x, Cl, Cf] = TrustRegion(StepBound(i), 0.1, 1e-9, 1000, ones(100,1), @cute_wrap, 0.5);
	err = norm(x-xout);
	display([err, Cl, Cf]);
end
