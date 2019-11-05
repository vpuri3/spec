function [z] = sem1dmesh(lx1,e,ifcheby)
%
% creates 1D sem mesh of poly. order lx1-1
% and e elements.
% z \in [0,1].
%

[z0,w0]=zwgll(lx1-1);
z0 = 0.5*(z0+1);      % [0,1]
w0 = 0.5*w0;

% get element mesh
if(ifcheby)
	dtheta=(0:e)'*pi/e;
	ze=(1-cos(dtheta))*0.5;
else
	ze=(0:1/e:1)'; % element mesh (size e+1)
end

z = kron(diff(ze),z0) + kron(ze(1:end-1),ones(lx1,1));
z = unique(z);

% w = kron(diff(ze),w0); % sum weights at element edges

%figure;plot(z,z*0,'kx');grid on;

end
