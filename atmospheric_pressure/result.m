xstart = [11, 3,0.5,0.025,−0.7,−1.3,0.04,−0.3,1.4]';

display('Midterm method')
tic; [x,k] = LeMa(100, 0.1, xstart, 1e-9, 1); toc
display('Cholesky method')
tic; [x,k] = LeMa(100, 0.1, xstart, 1e-9, 2); toc

display('Midterm method')
tic; [x1, k1, e1] = LeMaQN(100, 0.1, xstart, 1e-9, 1, x); toc
display('Cholesky method')
tic; [x2, k2, e2] = LeMaQN(100, 0.1, xstart, 1e-9, 2, x); toc

L1 = log(e1(2:end)./e1(1:end-1));
L2 = log(e2(2:end)./e2(1:end-1));

FNP1 = figure;
subplot(2,2,1);
plot(e1, 'k-o','MarkerSize', 4);
title('Eigenvalue decomposition method');
axis([0 40 -1 3]);
xlabel('Step number');
ylabel('Residual(i)');
subplot(2,2,2);
plot(L1, 'k-o','MarkerSize', 4);
title('Eigenvalue decomposition method');
axis([0 40 -5 1]);
xlabel('Step number');
ylabel('log(Residual(i+1) / Residual(i))');
subplot(2,2,3);
plot(e2, 'k-*','MarkerSize', 4);
title('Cholesky factorization method');
axis([0 40 -1 3]);
xlabel('Step number');
ylabel('Residual(i)');
subplot(2,2,4);
plot(L2, 'k-*','MarkerSize', 4);
title('Cholesky factorization method');
axis([0 40 -5 1]);
xlabel('Step number');
ylabel('log(Residual(i+1) / Residual(i))');
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'Figure1');
print(FNP1, '-dpdf', 'FNP1.pdf');


[ensoFit, ensoRes, ensoData] = ensoFitRes(x2);

FNP2 = figure;
subplot(2,2,1);
plot(ensoFit, ensoData, 'ks', 'MarkerSize', 3);
mmline = refline([1 0]);
set(mmline,'Color', 'r');
title('Fitted values versus actual data');
xlabel('Fitted value');
ylabel('Actual observation');
subplot(2,2,2);
plot(ensoFit, ensoRes, 'kd', 'MarkerSize', 3);
xlim([4 17])
title('Fitted values versus residuals');
xlabel('Fitted value');
ylabel('Residual');
subplot(2,2,3);
plot(ensoFit, abs(ensoRes), 'k*', 'MarkerSize', 3);
title('Fitted values versus absolute residuals');
xlim([3 17])
xlabel('Fitted value');
ylabel('Absolute residual');
subplot(2,2,4)
qqplot(ensoRes);
title('QQ plot of residuals versus standard normal distribution')
xlabel('Standard normal quantile');
ylabel('Sample quantile of residuals');
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'Figure2');
print(FNP2, '-dpdf', 'FNP2.pdf');
