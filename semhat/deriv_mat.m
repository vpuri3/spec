      function[d] =  deriv_mat1(x)
%
%     Compute the Lagrange interpolation matrix from x to x
%
      
      ni = length(x);
      a  = ones(ni,1);
      for i=1:ni;
        for j=1:(i-1);  a(i)=a(i)*(x(i)-x(j)); end;
        for j=(i+1):ni; a(i)=a(i)*(x(i)-x(j)); end;
      end;
      a=1./a; % These are the alpha_i's

      for j=1:ni; for i=1:ni; d(i,j)=x(i)-x(j); end; d(j,j)=1; end;
      d=1./d;
      for i=1:ni; d(i,i)=0; d(i,i)=sum(d(i,:)); end;

      for j=1:ni; for i=1:ni;
         if i~=j; d(i,j) = a(j)/( a(i)*(x(i)-x(j)));end;
      end;end;

