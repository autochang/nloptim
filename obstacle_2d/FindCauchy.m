function [xcauchy, countcauchy, breakpoint, gradient]= FindCauchy(gd,x,l,u,G,c)
	n=length(x); t=zeros(n,1); tbar=Inf(n,1); countcauchy=0;
	%gd = G*x+c; % steepest descent direction
	for i=1:n
		if gd(i)<0 & u(i)<Inf
			tbar(i) = (x(i)-u(i)) / gd(i);
		elseif gd(i)>0 & l(i)>-Inf
			tbar(i) = (x(i)-l(i)) / gd(i);
		end
	end
	tbartmp = tbar;
	tbar(tbar==0) = []; tt = [0; unique(tbar)]; % eliminate the duplicate and zero from tbar
	% compute the case t=0
	f0 = c'*x + x'*G*x/2;
	p0 = p1672(tbartmp,tt(1),gd);
	p0idx = find(p0 ~= 0);
	f1 = c(p0idx,1)'*p0(p0idx,1) + x'*G(:,p0idx)*p0(p0idx,1);
	f2 = p0'*G(:,p0idx)*p0(p0idx,1);
	% start searching cauchy point
	for j=1:length(tt)-1
		countcauchy = countcauchy + 1; 
		tDelta = -f1/f2;
		%display([f1 f2 tDelta]);
		x0 = plp(tt(j), x, gd, tbartmp);
		if f1>0
			xcauchy = x0; %=plp(tt(j), x, gd, tbar);
			breakpoint = tt(j); 
			gradient = f1;
			return
		elseif tDelta>=0 & tDelta<tt(j+1)-tt(j)
			xcauchy = plp(tt(j)+tDelta, x, gd, tbartmp);
			breakpoint = tt(j)+tDelta;
			gradient = f1;
			return
		end
		p1 = p1672(tbartmp, tt(j+1),gd);
		pDelta = p1-p0;
               	idx = find(pDelta ~= 0);
                comp0 = p0' * G(:,idx) * pDelta(idx, 1);
                comp1 = (c(idx,1)' + x0' * G(:,idx)) * pDelta(idx, 1) + (tt(j+1)-tt(j))*comp0;
                comp2 = 2*comp0 + pDelta' * G(:,idx) * pDelta(idx, 1);
                f0 = f0 + (tt(j+1)-tt(j))*f1 + (tt(j+1)-tt(j))^2*f2/2;
                f1 = f1 + (tt(j+1)-tt(j))*f2 + comp1;
                f2 = f2 + comp2;
		p0 = p1;
	end
	
