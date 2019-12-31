%
function [ux,uy] = grad(u,Dr,Ds,rx,ry,sx,sy);

	ur = ABu([],Dr,u);
	us = ABu(Ds,[],u);

	ux = rx .* ur + sx .* us;
	uy = ry .* ur + sy .* us;
end
