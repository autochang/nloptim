function midrad(c)
%MIDRAD       Display of interval gradients
%
%   midrad(c)
%

% written  10/16/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
%                                    take care of Matlab sparse Inf/NaN bug
% modified 02/11/06     S.M. Rump  SparseInfNanFlag removed
%

  global INTLAB_GRADIENT_NUMVAR

  loose = strcmp(get(0,'FormatSpacing'),'loose');

  numvar = size(c.dx,2);
  if numvar~=INTLAB_GRADIENT_NUMVAR
    warning('**** number of dependent variables and partial derivatives do not coincide')
  end

  sizecdx = size(c.x);
  if length(sizecdx)==2 & sizecdx(2)==1
    sizecdx = sizecdx(1);
  end

  if isa(c.x,'intval')
    if loose, disp(' '); end
    disp([ 'intval gradient value ' inputname(1) '.x = ' ])
    midrad(c.x,'',1)
    if loose, disp(' '); end

    if loose, disp(' '); end
    disp([ 'intval gradient derivative(s) ' inputname(1) '.dx = ' ])
    midrad( reshape( c.dx , [sizecdx INTLAB_GRADIENT_NUMVAR] ) , '' , 1 );
    if loose, disp(' '); end
  else
    display(c,inputname(1))
  end

  