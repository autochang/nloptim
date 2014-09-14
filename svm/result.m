% load the a9a dataset from csv file
filename = 'MatrixMat.csv';
MM = csvread(filename);

% draw the plots
ProblemSize = [10 20 40 80 100 200 400];
for i = 1:7
	n = ProblemSize(i);
	Data = MM(1:n,:);
	y = Data(:,1); H = MatQPD(Data);
	tic; [xsol, NumSV] = BCAugLag(ones(n,1), 10, 1e-4, 1e-4, 1000, Data); t(i)=toc;
	s(i)=NumSV;
end
save('SuVe.mat','t','s');

FNP4=figure;
subplot(1,2,1);
plot(ProblemSize./10, t, 'k-d');
xlim([0 41])
xlabel('n/10, n is the problem size');
ylabel('Run time (unit: second)');
subplot(1,2,2);
plot(ProblemSize./10, s./ProblemSize, 'k-s');
xlim([0 41])
xlabel('n/10, n is the problem size')
ylabel('Fraction of support vectors')
annotation('textbox', [0 0.9 0.08 0.04], 'String', 'Figure4');
print(FNP4, '-dpdf', 'FNP4.pdf');



