%
%	(v,b0*u -visc*\del^2 u)
%
function [Hu] =  hlmhltz(u,visc,b0,M,Qx,Qy,B,Dr,Ds,G11,G12,G22)

	Hu =      visc*lapl(u,[],[],[],Dr,Ds,G11,G12,G22);
	Hu = Hu +   b0*mass(u,B,[],[],[]);

	Hu = mass(Hu,[],M,Qx,Qy);
end
