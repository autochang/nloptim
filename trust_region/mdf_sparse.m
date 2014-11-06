function [F] = mdf_sparse(Q,lambda,delta)
tau = (delta-lambda) .* (lambda<delta);
Q1 = sparse(Q); 
F = Q1 * sparse(diag(tau)) * Q1';
