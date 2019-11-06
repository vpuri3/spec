% Real fft
%
%
%   Given u(1:n), this returns the vector uh(1:n).
%
%   If
%             m_c                              m_s
%      u(i) = sum  a_k cos( 2 pi k (i-1)/n ) + sum  b_k sin( 2 pi k (i-1)/n )
%             k=0                              k=0
%
%   then the coefficients in uh() are given by
%
%   uh:
%
%     a0   a1     a2     a3     a4     ...     a_m_c
%             b1     b2     b3     b4     ...     b_m_s
%
%   where
%
%           m_s = m_c = n/2 if n is odd
%
%           m_s = m_c-1, m_c = n/2 if n is even
%
%   To invert and recover u(), use  u=irfft(uh);
%
%   If u is an matrix (n x m) matrix, this routine returns the
%   transformed columns uh(:,j), j=1:m.
%
%

function [uh] = rfft(u);

n  = size(u,1);
m  = size(u,2);
n2 = floor(n/2);
n1 = n2-1;
if mod(n,2)==1; n1=n2; end;

a = zeros(n,m);
b = zeros(n,m);
z = fft(u);

a0 = z(1,:)/n;
a(1:n-1,:) =  2*real(z(2:end,:))/n;
b(1:n-1,:) = -2*imag(z(2:end,:))/n;

r=zeros(n,m);
uh(1,:) = a0;
uh(2:2:n-1,:) =  a(1:n1,:);
uh(3:2:n  ,:) =  b(1:n1,:);
if mod(n,2)==0, uh(n,:) = a(n2,:)/2; end;


