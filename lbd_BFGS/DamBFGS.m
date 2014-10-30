function [bs, theta, yhat] = DamBFGS(S, Bs, Theta, Yhat, s, y, B0)
	[n, m]=size(S); bs=0;
	for j=1:m
		C1 = dot(Bs(:,j),s) / dot(S(:,j),Bs(:,j));
		C2 = dot(Yhat(:,j),s) / dot(S(:,j),Yhat(:,j));
		bs = bs - C1*Bs(:,j) + C2*Yhat(:,j);
	end
	bs = bs + B0*s;
	if dot(s,y) >= 0.2*dot(s, bs)
		theta =1;
	else
		theta = 0.8*dot(s, bs) / (dot(s, bs)-dot(s, y));
	end
	yhat = theta*y + (1-theta)*bs;
