function [z,w] = semmesh(E,lx1,ifcheby)
%
% creates 1D sem mesh of poly. order lx1-1
% and e elements [-1,1].
%

[z0,w0]=zwgll(lx1-1);
z0=0.5*(z0+1); % [0,1]
w0=0.5*w0;

% element mesh
if(ifcheby)
	dtheta = (0:E)'*pi/E;
	ze = -cos(dtheta);
else
	ze=(-1:2/E:1)'; % element mesh (size e+1)
end

z = kron(diff(ze),z0) + kron(ze(1:end-1),ones(lx1,1));
%z = unique(z);

w = kron(diff(ze),w0);
%Q = semq(E,lx1-1,0);
%w = Q'*w;

%figure;plot(z,z*0,'kx');grid on; sum(w)

end
