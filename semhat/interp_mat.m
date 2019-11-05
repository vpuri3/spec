      function[J] =  interp_mat(xo,xi)
%
%     Compute the Lagrange interpolation matrix from xi to xo
%
      
      no = length(xo);
      ni = length(xi);
      a  = ones(ni,1);
      for i=1:ni;
        for j=1:(i-1);  a(i)=a(i)*(xi(i)-xi(j)); end;
        for j=(i+1):ni; a(i)=a(i)*(xi(i)-xi(j)); end;
      end;
      a=1./a;

      J = zeros(no,ni);
      s = ones(ni,1); t = ones(ni,1);
      for i=1:no;
        x=xo(i);
        for j=2:ni; s(j)      = s(j-1)*(x-xi(j-1));
                    t(ni+1-j) = t(ni+2-j)*(x-xi(ni+2-j)); end;
        J(i,:)=a.*s.*t;
      end;
