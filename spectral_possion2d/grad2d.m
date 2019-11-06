%
function [ux,uy] = grad2d(u,Ir,Is,Dr,Ds,rx,ry,sx,sy);
	
	ur = ABu(Is,Dr,u);
	us = ABu(Ds,Ir,u);

	ux = ur.*rx + us.*sx;
	uy = ur.*ry + us.*sy;
