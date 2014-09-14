FNP3=figure;
NN = [16 24 32 40];
for i=1:4
	N = NN(i);
	[G,l,f,X,Y]=obstacle_mihai(N);
	c = -f;
	[xsol, k1, k2, k3] = GPCG(l, G, c, l, Inf(N*N,1), 2, 3, 1e-9);
	subplot(2,2,i);
	displayObst(N,X,Y,xsol);
	title(['N=' num2str(N)]);
end
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'Figure3');
print(FNP3, '-dpdf', 'FNP3.pdf');

