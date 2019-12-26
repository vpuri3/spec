%
%	(v,b0*u -visc*\del^2 u)
%
function [Hu] =  hlmhltz(u,visc,b0,B,Ir,Is,Dr,Ds,G11,G12,G22)

	Hu =      visc*lapl(u,Ir,Is,Dr,Ds,G11,G12,G22);
	Hu = Hu +   b0*mass(u,B,Ir,Is);

end
