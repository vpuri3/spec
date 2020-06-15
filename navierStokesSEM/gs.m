%
% Gather-Scatter
% 
function [Gu] = gs(u,Qx,Qy);

	Qu = ABu(Qy',Qx', u); % gather
	Gu = ABu(Qy ,Qx ,Qu); % scatter

end
