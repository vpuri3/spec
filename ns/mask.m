%
% masks dirichlet boundary points
%
function [v] = mask(M,u)

	%v = ABu(Ry'*Ry,Rx'*Rx,u);
	v = M .* u;

end
