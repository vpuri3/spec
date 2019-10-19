%
function [ux,uy] = gradu(I,D,rx,ry,sx,sy,u);

	ux = rx.*ABu(I,D,u) + sx.*ABu(D,I,u);
	uy = ry.*ABu(I,D,u) + sy.*ABu(D,I,u);
