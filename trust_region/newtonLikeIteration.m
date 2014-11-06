function [x,gradNorm]=newtonLikeIteration(functionHandle,xStartPoint,iterationType,iterIndex)
% computes one iteration of Newton Like  method with a diagonal
% perturbation
% INPUT:    functionHandle:     (pointer) Function defining problem
%           xStartPoint:        (vector)  The Starting Point
%           iterIndex:          (integer) The index of the iterate
%           iterationType:      k=1:    Newton's Method
%                               k=2:    Hessian peturbed by identity, I
%                               k=3:    Hessian peturbed by o(iterInd)*I
% OUTPUT:   x:                  (vector) Next iteration Point. 


% test for consistency
if iterIndex <=0
    % the iteration cannot be at $0$, throw an error
    error('The iteration index cannot be 0');
end

% initialize
x=xStartPoint;

% decide which perturbation you want to use

switch iterationType
    case 1,
        diagonalPerturbation=0;
    case 2, 
        diagonalPerturbation=1;
    case 3, 
        diagonalPerturbation=1/iterIndex;
    otherwise,
        error('That iteration option is not implemented');
end

% get the function data
[f,g,H]=functionHandle(x,2);

% the actual calculation
gradNorm=norm(g);
H=H+diagonalPerturbation*eye(length(x));
x=x-H\g;


        



