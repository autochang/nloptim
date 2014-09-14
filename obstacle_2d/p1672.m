function p = p1672(tbar,t,gd)
	n=length(tbar); p=zeros(n,1);
	for i=1:n
		if tbar(i)>t
			p(i)=-gd(i);
		end	
	end
