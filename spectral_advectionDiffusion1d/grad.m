%
function [ux] = grad(u,Dx);

	ux = Dx * ux;
