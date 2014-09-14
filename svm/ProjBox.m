function Proj = ProjBox(g, l, u)
	n = length(g); Proj=zeros(n,1);
	for i=1:n
		if g(i)<=l(i)
			Proj(i)=l(i);
		elseif g(i)>=u(i)
			Proj(i)=u(i);
		else
			Proj(i)=g(i);
		end
	end
