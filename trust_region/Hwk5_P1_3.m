xout = newtonLikeMethod(@cute_wrap, ones(100,1), 1, 1e-9, 1000);
[xT] = TrustRegionSol(5, 0.1, 1000, ones(100,1), @cute_wrap, 0.5);
[xL] = LineSearchSol(1, 0.5, 0.5, 1000, ones(100,1), @cute_wrap, 0.5);

[EQT] = SuperLinear(xT(:, 1:20), xT(:,1000));
[EQL] = SuperLinear(xL(:, 1:20), xL(:,1000));

((2:100).^(-(2:100))) ./ ((1:99).^(-(1:99)))


