clear;
run('Intlab_V5.5/startup.m');

ProblemSize = [100 300 1000 3000];
M = [10 20 40 80 100];

K = zeros(1,5); Time = zeros(1,5);	
for I = 1:5	
	tic; [x, K(I)] = LBD_BFGS(-ones(ProblemSize(I),1), 100, 1e-9, @Rosenbrock_wrap); Time(I)=toc;
end

NumIterAVG = zeros(5,4); NumIterSTD = zeros(5,4);
RunTimeAVG = zeros(5,4); RunTimeSTD = zeros(5,4);

for I = 1:4 % index problem size
	for J = 1:5 % index storage pair
		KK = zeros(10,1); TT = zeros(10,1);
		for L = 1:10
			N = ProblemSize(I);
			tic; [x, KK(L)] = LBD_BFGS(-ones(N,1)+randn(N,1)*0.05, M(J), 1e-4, @Rosenbrock_wrap); TT(L)=toc;
		end
		NumIterAVG(J,I) = mean(KK); NumIterSTD(J,I)=std(KK);
		RunTimeAVG(J,I) = mean(TT); RunTimeSTD(J,I)=std(TT);
	end
end
NumIter(:,1:2:7)=NumIterAVG; NumIter(:,2:2:8)=NumIterSTD;
RunTime(:,1:2:7)=RunTimeAVG; RunTime(:,2:2:8)=RunTimeSTD;
save('H7P2b.mat','NumIterAVG','NumIterSTD','RunTimeAVG', 'RunTimeSTD');
display(NumIter)
display(RunTime)

KKK = zeros(1,5); TTT = zeros(1,5);
TE = 10.^[-2 -3 -4];      
for I = 1:5
        tic; [x, KKK(I)] = LBD_BFGS(-ones(3000,1), 100, TE(I), @Rosenbrock_wrap); TTT(I)=toc;
end
display(KKK)
display(TTT)


NumIterAVG = zeros(5,3); NumIterSTD = zeros(5,3);
RunTimeAVG = zeros(5,3); RunTimeSTD = zeros(5,3);

for I = 1:3 % index error
        for J = 1:5 % index storage pair
                KK = zeros(10,1); TT = zeros(10,1);
                for L = 1:10
                         tic; [x, KK(L)] = LBD_BFGS(-ones(300,1)+randn(300,1)*0.05, M(J), TE(I), @Rosenbrock_wrap); TT(L)=toc;
                end
                NumIterAVG(J,I) = mean(KK); NumIterSTD(J,I)=std(KK);
                RunTimeAVG(J,I) = mean(TT); RunTimeSTD(J,I)=std(TT);
        end
end
NumIter(:,1:2:5)=NumIterAVG; NumIter(:,2:2:6)=NumIterSTD;
RunTime(:,1:2:5)=RunTimeAVG; RunTime(:,2:2:6)=RunTimeSTD;
save('H7P2c.mat','NumIterAVG','NumIterSTD','RunTimeAVG', 'RunTimeSTD');
display(NumIter)
display(RunTime)
