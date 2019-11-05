      function [z,w] = zwglc(p);

%
%
% computes the p+1 Gauss-Lobatto-Chebyshev nodes z on [-1,1]
% i.e. the zeros of the first derivative of the Legendre polynomial
% of degree p plus -1 and 1
%
% and the p+1 weights w
%


n = p+1;

z(1:n)=0;
w(1:n)=0;

dt = -pi/p;
theta = pi:dt:0; theta=theta';
z=cos(theta);
w=(1/p) + 0*z; 
w(1)=0.5*w(1);
w(n)=0.5*w(n);


