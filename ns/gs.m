%
% Gather-Scatter
% 
function [Gu] = gs(u,Qx,Qy);

if(length(Qx)==0 & length(Qy)==0);
	Gu = u;
else
	Qu = ABu(Qy',Qx', u); % gather
	Gu = ABu(Qy ,Qx ,Qu); % scatter
end

end
