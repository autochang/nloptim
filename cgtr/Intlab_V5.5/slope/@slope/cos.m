function u = cos(a)
%COS          Slope cosine cos(a)
%

% written  12/06/98     S.M. Rump
% modified 04/04/04     S.M. Rump  set round to nearest for safety
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  global INTLAB_SLOPE

  u = a;

  u.r = cos(a.r);
  indexc = 1:INTLAB_SLOPE.NUMVAR;
  indexr = 2:INTLAB_SLOPE.NUMVAR+1;
  Xxs = hull(a.r(:,indexc),a.r(:,indexr));
  u.s = - sin(Xxs) .* a.s;
  
  if rndold~=0
    setround(rndold)
  end
