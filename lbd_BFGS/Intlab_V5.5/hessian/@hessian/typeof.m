function t = typeof(a)
%TYPEOF       Type of a
%
%   t = typeof(a)
%
%For details, see intval\typeof and intval\typeadj.
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

% hessian\@hessian\typeof:  input  a  must be hessian
  if isa(a.x,'intval')
    t = 'hessianintval';
  else
    t = 'hessian';
  end
