function [x] = LineSearchSol(alpha0, c, rho, NumIter, xstart, FHandle, delta)
n = length(xstart);
pd = zeros(n,1);
x = zeros(n,NumIter);
x(:,1) = xstart;  
for k = 1:(NumIter-1)
        [f, g, H] = FHandle(x(:,k), 2);
        %[L, D, P] = ldl(H);
        [L, D, p] = ldl(H, 'vector'); % ldl
        [Q, S, Q] = svd(full(D)); % better than svds, which is better for a few svs
	%[Q, S, Q] = svds(D, n);% evd
        F = mdf_sparse(Q, diag(S), delta);
        %B = P'*L*(D+F)*L'*P; % modified Hessian
        %pd = -B \ g; % modified Newton direction
        pd(p, :) = -L' \ ((D+F) \ (L \ g(p, :) ));
        [alpha, fct]= bkt(x(:,k), pd, FHandle, alpha0, c, rho); % backtracking
	x(:,k+1) = x(:,k) + alpha * pd;
end

