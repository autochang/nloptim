function L = NextPowerTwo(p)
%NEXTPOWERTWO Smallest power of 2 not less than abs(p) for p~=0
%
%   L = NextPowerTwo(p)
%
%On return, abs(p) <= 2^L and L is smallest possible.
%
%Implements Algorithm 3.5 from
%  S.M. Rump, T. Ogita, S. Oishi: Accurate Floating-point Summation I: 
%    Faithful Rounding, to appear. 
%Requires 4 flops w/o scaling.
%
%Reference implementation! Slow due to interpretation!
%
%Improvement: in the following simpler version, the "if"-statement is omitted
%

% written  03/03/07     S.M. Rump
% modified 09/19/07     S.M. Rump  improved implementation: omit branch
%

% Improved algorithm w/o branch
  feature accel off                 % switch off Matlab code optimization
  if isa(p,'double'), prec='double'; else prec='single'; end
  u = eps(prec)/2;
  if abs(p)<pred(u*realmax(prec))   % no large input
    q = u^(-1) * p;
    L = abs( ( q-p ) - q );
  else                              % large input
    L = 2*u^(-1)*NextPowerTwo(0.5*u*p);
  end
  
return

% The algorithm from the cited paper
%   feature accel off                 % switch off Matlab code optimization
%   if isa(p,'double'), prec='double'; else prec='single'; end
%   u = eps(prec)/2;
%   if abs(p)<pred(u*realmax(prec))   % no large input
%     q = u^(-1) * p;
%     L = abs( ( q+p ) - q );
%     if L==0
%       L = abs(p);
%     end
%   else                              % large input
%     L = 2*u^(-1)*NextPowerTwo(0.5*u*p);
%   end
