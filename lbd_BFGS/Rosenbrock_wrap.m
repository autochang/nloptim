function [f,g,H]=Rosenbrock_wrap(x,nrderivs);
% INPUT:    x:          (vector) input value 
%           nrderivs:   (integer) the number of derivs needed. 
%                       0: return only function 
%                       1: return function  and gradient.
%                       2+: return function, gradient, and Hessian H. 
% OUTPUT:   f:          (scalar) Function value
%           g:          (vector) Gradient
%           H:          (matrix) Hessian

% consistency tests

    
% choose actual derivative info needed
switch nrderivs
    case 0,
        f=Rosenbrock(x);
    case 1,
        gradientBundle=Rosenbrock(gradientinit(x));
        f=gradientBundle.x;
        g=gradientBundle.dx';
    case 2, 
        hessianBundle=Rosenbrock(hessianinit(x)); 
        f=hessianBundle.x;
        g=hessianBundle.dx';
        H=hessianBundle.hx;
    otherwise,
        % throw an error, case not defined
        error('This derivative option is not defined');
end

        
