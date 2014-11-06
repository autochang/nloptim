delta = 10.^[-2;-1;0;1;2;3];
xout = newtonLikeMethod(@cute_wrap, ones(100,1), 1, 1e-9, 1000);
for i = 1:6
	[x, Cl, Cf] = TrustRegion(1, 0.1, 1e-9, 1000, ones(100,1), @cute_wrap, delta(i));
	err = norm(x-xout);
	display([err, Cl, Cf]);
end
