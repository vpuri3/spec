%
%     (v,-\del^2 u)
%
function [Au] = lapl(u,M,Qx,Qy,Dr,Ds,G11,G12,G22)

uu = mask(u,M);
uu = ABu(Qy,Qx,uu);

ur = ABu([],Dr,uu);
us = ABu(Ds,[],uu);

wr = G11.*ur + G12.*us;
ws = G12.*ur + G22.*us;

Au = ABu([],Dr',wr) + ABu(Ds',[],ws);

Au = ABu(Qy',Qx',Au);
Au = mask(Au,M);

end
