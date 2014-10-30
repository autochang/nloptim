clear;
run('Intlab_V5.5/startup.m');

ProbSize = [100 300 1000 3000];
NumIter = zeros(1,4); NumIterCG = zeros(1,4);
LDLNumIter = zeros(1,4); LDLNumIterCG = zeros(1,4);
RunTime = zeros(1,4); LDLRunTime = zeros(1,4);

H8P1 = figure;
for i = 1:4
	n = ProbSize(i);
	t0=tic; [x, NumIter(i), NumIterCG(i), e0]=TR_HWK8(100, 0.1, 1e-9, randn(n,1), @Rosenbrock_wrap); RunTime(i)=toc(t0);
	t1=tic; [x, LDLNumIter(i),LDLNumIterCG(i), e1]=LDLTR_HWK8(100, 0.1, 1e-9,randn(n,1), @Rosenbrock_wrap); LDLRunTime(i)=toc(t1);
	subplot(2, 2, i);
	plot(e0, 'k-o','MarkerSize', 3.5);
	hold on;
	plot(e1, 'k-*','MarkerSize', 3.5);
	hold off;
	PT = strcat('n=',num2str(ProbSize(i)));
	title(PT);
	xlabel('Step number');
	ylabel('Residual')
	mleg = legend('Without preconditioning', 'LDL preconditioning');
	set(mleg, 'Location', 'northeast');
end
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'FigureR1');
print(H8P1, '-dpdf', 'H8P1R.pdf');

H8P2 = figure;
for i = 1:4
        n = ProbSize(i);
        t0=tic; [x, NumIter(i), NumIterCG(i), e0]=TR_HWK8(100, 0.1, 1e-9, randn(n,1), @Rosenbrock_wrap); RunTime(i)=toc(t0);
        t1=tic; [x, LDLNumIter(i),LDLNumIterCG(i), e1]=LDLTR_HWK8(100, 0.1, 1e-9,randn(n,1), @Rosenbrock_wrap); LDLRunTime(i)=toc(t1);
        subplot(2, 2, i);
        plot(e0(2:end)./e0(1:end-1), 'k-o','MarkerSize', 3.5);
        hold on;
        plot(e1(2:end)./e1(1:end-1), 'k-*','MarkerSize', 3.5);
        hold off;
        PT = strcat('n=',num2str(ProbSize(i)));
        title(PT);
 	axis([0 20 0 1.1]);
        xlabel('Step number');
        ylabel('Residual(i+1) / Residual(i)')
        mleg = legend('Without', 'LDL');
        set(mleg, 'Location', 'southwest');
end
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'FigureR2');
print(H8P2, '-dpdf', 'H8P2R.pdf');


H8P3 = figure;
for i = 1:4
        n = ProbSize(i);
        t0=tic; [x, NumIter(i), NumIterCG(i), e0]=TR_HWK8(100, 0.1, 1e-9, randn(n,1), @Rosenbrock_wrap); RunTime(i)=toc(t0);
        t1=tic; [x, LDLNumIter(i),LDLNumIterCG(i), e1]=LDLTR_HWK8(100, 0.1, 1e-9,randn(n,1), @Rosenbrock_wrap); LDLRunTime(i)=toc(t1);
        subplot(2, 2, i);
        plot(log(e0(2:end)./e0(1:end-1)), 'k-o','MarkerSize', 3.5);
        hold on;
        plot(log(e1(2:end)./e1(1:end-1)), 'k-*','MarkerSize', 3.5);
        hold off;
        PT = strcat('n=',num2str(ProbSize(i)));
        title(PT);
	axis([0 20 -20 1]);
        xlabel('Step number');
        ylabel('log( Residual(i+1)/Residual(i) )')
        mleg = legend('Without', 'LDL');
        set(mleg, 'Location', 'southwest');
end
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'FigureR3');
print(H8P3, '-dpdf', 'H8P3R.pdf');
