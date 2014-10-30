clear;
run('Intlab_V5.5/startup.m');

ProbSize = [100 300 1000 3000];
NumIter = zeros(1,4); NumIterCG = zeros(1,4);
SDNumIter = zeros(1,4); SDNumIterCG = zeros(1,4);
LDLNumIter = zeros(1,4); LDLNumIterCG = zeros(1,4);
SDLDLNumIter = zeros(1,4); SDLDLNumIterCG = zeros(1,4);
RunTime = zeros(1,4); LDLRunTime = zeros(1,4);
SDRunTime = zeros(1,4); SDLDLRunTime = zeros(1,4);

for i = 1:4
        n = ProbSize(i);
	k = zeros(1,20); kk = zeros(1,20);
	LDLk = zeros(1,20); LDLkk = zeros(1,20);
	t = zeros(1,20); LDLt = zeros(1,20);
	for j = 1:20
        	t0=tic; 
		[x, k(j), kk(j), e0]=TR_HWK8(100, 0.1, 1e-9, randn(n,1), @Rosenbrock_wrap); 
		t(j)=toc(t0);
        	t1=tic; 
		[x, LDLk(j),LDLkk(j), e1]=LDLTR_HWK8(100, 0.1, 1e-9, randn(n,1), @Rosenbrock_wrap); 
		LDLt(j)=toc(t1);
	end
	NumIter(i)=mean(k); NumIterCG(i)=mean(kk);
	SDNumIter(i)=std(k); SDNumIterCG(i)=std(kk);
	LDLNumIter(i)=mean(LDLk); LDLNumIterCG(i)=mean(LDLkk);
	SDLDLNumIter(i)=std(LDLk); SDLDLNumIterCG(i)=std(LDLkk);
	RunTime(i)=mean(t); LDLRunTime(i)=mean(LDLt);
	SDRunTime(i)=std(t); SDLDLRunTime(i)=std(LDLt);
end

TIME(1:2:7)=RunTime; TIME(2:2:8)=SDRunTime;
LDLTIME(1:2:7)=LDLRunTime; LDLTIME(2:2:8)=SDLDLRunTime;

K(1:2:7)=NumIter; K(2:2:8)=SDNumIter;
LDLK(1:2:7)=LDLNumIter; LDLK(2:2:8)=SDLDLNumIter;

KK(1:2:7)=NumIterCG; KK(2:2:8)=SDNumIterCG;
LDLKK(1:2:7)=LDLNumIterCG; LDLKK(2:2:8)=SDLDLNumIterCG;

display(TIME)
display(LDLTIME)
display(K)
display(LDLK)
display(KK)
display(LDLKK)
