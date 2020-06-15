%
% (v,b1*u - visc*u_xx)
%
function [Hu] =  hmhltz(u,b0,visc,B,D)
	Hu =        b0*mass(u,B);
	Hu = Hu + visc*lapl(u,B,D);
end
