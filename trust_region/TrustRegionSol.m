function [x] = TrustRegionSol(StepBound, eta, NumIter, xstart, FHandle, delta)
n = length(xstart);
Delta = StepBound / 2;
pd = zeros(n,1);
x = zeros(n,NumIter);
x(:,1) = xstart;
for k = 1:(NumIter-1)
        [f, g, H] = FHandle(x(:, k), 2); % one function eval
        [L, D, p] = ldl(H, 'vector');
        [Q, S, Q] = svd(full(D));
        F = mdf_sparse(Q, diag(S), delta);
        B(p,p) = L * (D+F) * L';
        pd = dogleg(B, g, Delta); % one linear solve
        rho = reductratio(x(:, k), pd, f, g, B, FHandle); % one function eval
        if rho < 0.25
                Delta = Delta / 4;
        elseif (rho > 0.75) & (norm(p) == Delta)
                Delta = min(2*Delta, StepBound);
        else
                Delta = Delta;
        end
        if rho > eta
                x(:,k+1) = x(:,k) + pd;
        else
                x(:,k+1) = x(:,k);
        end
end

