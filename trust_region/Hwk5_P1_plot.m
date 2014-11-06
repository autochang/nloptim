StopError = 1e-9; MaxIter = 1000;

TRCl = zeros(20,1); TRCf = zeros(20,1);
LSCl = zeros(20,1); LSCf = zeros(20,1);

for j = 1:20
	n = 100 * j;
	[xTR, TRCl(j), TRCf(j)] = TrustRegion(5, 0.1, 1e-9, 1000, ones(n,1), @cute_wrap, 0.5);
	[xLS, LSCl(j), LSCf(j)] = LineSearch(1, 0.5, 0.5, 1e-9, 1000, ones(n,1), @cute_wrap, 0.5);
end

H = figure;

subplot(2, 1, 1);
index = 100*(1:20);
plot(index, TRCl, 'k-o', index, LSCl, 'k-*');
xlim([0 2100]);
title('Number of linear systems solved');
xlabel('Problem size');
mleg1 = legend('Dogleg trust region', 'Backtracking line search');
set(mleg1, 'Location', 'northeast');

subplot(2,1,2);
index = 100*(1:20);
plot(index, TRCf, 'k-o', index, LSCf, 'k-*');
xlim([0 2100]);
title('Number of functions evaluated');
xlabel('Problem size');
mleg2 = legend('Dogleg trust region', 'Backtracking line search');
set(mleg2, 'Location', 'northeast');
 
print(H,'-dpdf', 'Hwk5_P1_plot.pdf');
