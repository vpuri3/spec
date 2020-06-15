%
% masks dirichlet boundary points
%
function [Mu] = mask(u,M)

	if(length(M)==0); Mu = u;
	else              Mu = M .* u;	%Mu = ABu(Ry'*Ry,Rx'*Rx,u);
	end

end
