function [xout,iteratesGradNorms]=newtonLikeMethod(functionHandle,xStartPoint,iterationType,stopTolerance,maxIterations)
%  [xout,iteratesGradNorms]=newtonLikeMethod(functionHandle,xStartPoint,iterationType,stopTolerance,maxIterations)
% PURPOSE:  computes the outcome of a Newton-like Method with a perturbed diagonal
% INPUT:    functionHandle:     (pointer) Handle to the optimization 
%                               problem to be solved
%           xStartPoint:        (vector) The starting point
%           iterationType:      (integer) The type of the diagonal peturbation to be
%                               used
%           stopTolerance:      (scalar) the gradient size at which the iteration
%                               will stop.
%           maxIterations:      (integer) The maximum number of iterations for which
%                               the algorithm should be run
% OUTPUT:   xout:               (vector) The final output
%           iteratesGradNorms:  (vector) The norms of the gradients

% consistency tests
if stopTolerance <=0
    error('Tolerance Parameter is invalid');
end

if maxIterations <=0
    error('The Number of Iterations is Invalid');
end

% initialization
iterCounter=1;
x=xStartPoint;
currentGradNorm=2*stopTolerance;

% main loop

while (currentGradNorm >= stopTolerance) & (iterCounter < maxIterations)
    [x,currentGradNorm]=newtonLikeIteration(functionHandle,x,iterationType,iterCounter);
    iteratesGradNorms(iterCounter,1)=currentGradNorm;
    iterCounter=iterCounter+1;
end

xout=x;
