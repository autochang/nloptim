function pd = LimBFGS(g, S, Y, H)
	q = g; 
	m = size(S,2); alpha = zeros(m, 1);
	for i = m:-1:1
		alpha(i) = dot(S(:,i), q) / dot(Y(:,i), S(:,i));
		q = q - alpha(i)*Y(:,i);	
	end
	r = H*q;
	for i = 1:m
		beta = dot(Y(:,i), r) / dot(Y(:,i), S(:,i));
		r = r + (alpha(i)-beta)*S(:,i);
	end
	pd = -r;
