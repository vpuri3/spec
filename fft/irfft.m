
%
% Inverse Real fft
%
%  See rfft.m for comments
%
%

function [u] = irfft(uh);

n  = size(uh,1);  m  = size(uh,2);
l  = n+1;
n2 = floor((n-1)/2);

u=zeros(n,m);
z=zeros(n,m);

z(1,:) = n*complex(uh(1,:),0);

for k=1:n2;
    z(1+k,:) = n*complex(uh(2*k,:),-uh(2*k+1,:))/2;
    z(l-k,:) = n*complex(uh(2*k,:), uh(2*k+1,:))/2;
end;
if mod(n,2)==0, z(n2+2,:) = n*complex(uh(n,:),0); end;

u = real( ifft(z) );
