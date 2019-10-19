      function[Bh,Dh,z,w] =  semhat(N)
%
%                                                 ^
%     Compute the single element 1D SEM Stiffness Mass, and Convection
%     matrices, as well as the points and weights for a polynomial
%     of degree N
%
      [z,w] = zwgll(N);

      Bh    = diag(w); % = l(x) * l(x)'
      Dh    = dhat(z); % l'(x) = Dh * l(x)

     %Ah    = Dh'*Bh*Dh; % = l'(x)' * l'(x)
     %Ch    = Bh*Dh;

