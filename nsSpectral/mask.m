%
% masks dirichlet boundary points
%
function [v] = mask(u,Rx,Ry)

	v = ABu(Ry'*Ry,Rx'*Rx,u);
	%v = msk .* u;

end
