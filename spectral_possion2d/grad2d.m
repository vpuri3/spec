%
function [ux,uy] = grad2d(u,I,D,rx,ry,sx,sy);
	
	ur = ABu(I,D,u);
	us = ABu(D,I,u);

	ux = ur.*rx + us.*sx;
	uy = ur.*ry + us.*sy;
