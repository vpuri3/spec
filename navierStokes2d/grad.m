%
function [ux,uy] = grad(u,Ir,Is,Dr,Ds,rx,ry,sx,sy);

	ur = ABu(Is,Dr,u);
	us = ABu(Ds,Ir,u);

	ux = rx .* ur + sx .* us;
	uy = ry .* ur + sy .* us;
end
