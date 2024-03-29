function r = sqrt(a)
%SQRT         Hessian (elementwise) square root
%

% written  04/04/04     S.M. Rump
% modified 04/06/05     S.M. Rump  rounding unchanged
%

  e = 1e-30;
  if 1+e==1-e                           % fast check for rounding to nearest
    rndold = 0;
  else
    rndold = getround;
    setround(0)
  end

  K = prod(size(a.x));
  if K==1                   % scalar hessian
    
    r.x = sqrt(a.x);
    rx1 = 0.5 / r.x;
    r.dx = a.dx * rx1;
    r.hx = ( a.hx - reshape( ((0.25/a.x)*a.dx) * a.dx.',size(a.hx)) ) * rx1;
    
  else                      % matrix hessian
    
    global INTLAB_HESSIAN_NUMVAR
    N = INTLAB_HESSIAN_NUMVAR;        
    N2 = N^2;
    
    r.x = sqrt(a.x);
    if issparse(a.hx)               % input sparse
      
      ax = 0.5./full(r.x(:));
      sizeax = length(ax);
      [ia,ja,sa] = find(a.dx);
      % check for emptyness: cures Matlab bug
      % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
      if isempty(ia)
        r.dx = sparse([],[],[],N,sizeax);
        r.hx = sparse([],[],[],N2,sizeax);
      else
        rdx = ax(ja).*sa(:);
        adx1 = (-0.25) ./ a.x(ja);
        if isa(a.x,'intval')          % sparse intval
          adx1 = times(adx1(:),rdx,0);
          if isreal(rdx)
            r.dx = intval( sparse(ia,ja,rdx.inf,N,sizeax) , sparse(ia,ja,rdx.sup,N,sizeax) , 'infsup' );
          else
            r.dx = intval( sparse(ia,ja,rdx.mid,N,sizeax) , sparse(ia,ja,rdx.rad,N,sizeax) , 'midrad' );
          end
          if adx1.complex
            adx1 = intval( sparse(ia,ja,adx1.mid,N,sizeax) , sparse(ia,ja,adx1.rad,N,sizeax) , 'midrad' );
          else
            adx1 = intval( sparse(ia,ja,adx1.inf,N,sizeax) , sparse(ia,ja,adx1.sup,N,sizeax) , 'infsup' );
          end
        else                          % sparse point  
          r.dx = sparse(ia,ja,rdx,N,sizeax);        
          adx1 = sparse(ia,ja,adx1(:).*rdx,N,sizeax);        
        end                           
        r.hx = adx2rhx(N,sizeax,adx1,a.dx);
      end
      [ia,ja,sa] = find(a.hx);      % sparse point or intval
      % check for emptyness: cures Matlab bug
      % a=sparse([],[],[],2,1), [i,j,s]=find(a), s(i).*s(:)  yields error
      if ~isempty(ia)
        if isa(a.x,'intval')
          rhx = times(ax(ja),sa(:),0);
          if rhx.complex
            r.hx = r.hx + intval( sparse(ia,ja,rhx.mid,N2,sizeax) , sparse(ia,ja,rhx.rad,N2,sizeax) , 'midrad' );
          else
            r.hx = r.hx + intval( sparse(ia,ja,rhx.inf,N2,sizeax) , sparse(ia,ja,rhx.sup,N2,sizeax) , 'infsup' );
          end
        else
          r.hx = r.hx + sparse(ia,ja,ax(ja).*sa(:),N2,sizeax);
        end
      end
      
    else                            % input full
      
      r.x = sqrt(a.x);
      rx1 = 0.5 ./ ( r.x(:).' );
      rx1 = rx1(ones(N*N,1),:);
      r.dx = a.dx .* rx1(1:N,:);
      adx = repmat(0.25./(a.x(:).'),N,1) .* a.dx;
      r.hx = ( a.hx - adx(repmat(1:N,N,1),:) .* a.dx(repmat(1:N,1,N),:) ) .* rx1;
      
    end
    
  end
  
  r = class(r,'hessian');
  
  if rndold~=0
    setround(rndold)
  end
