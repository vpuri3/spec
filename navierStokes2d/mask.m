%----------------------------------------------------------------------
function [v] = mask(u,msk)
	%v = ABu(Ry'*Ry,Rx'*Rx,u);
	v = msk .* u;
end