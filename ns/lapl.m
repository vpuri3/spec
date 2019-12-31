%
%     (v,-\del^2 u)
%
function [Au] = lapl(u,Qx,Qy,Dr,Ds,G11,G12,G22)

ul = ABu(Qy,Qx,u);

ur = ABu([],Dr,ul);
us = ABu(Ds,[],ul);

wr = G11.*ur + G12.*us;
ws = G12.*ur + G22.*us;

Au = ABu([],Dr',wr) + ABu(Ds',[],ws);

Au = ABu(Qy',Qx',Au);

end

