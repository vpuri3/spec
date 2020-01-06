%
%     (v,-\del^2 u)
%
function [Au] = lapl(u,M,Qx,Qy,Dr,Ds,G11,G12,G22)

uu = mask(u,M);

ur = ABu([],Dr,uu);
us = ABu(Ds,[],uu);

wr = G11.*ur + G12.*us;
ws = G12.*ur + G22.*us;

Au = ABu([],Dr',wr) + ABu(Ds',[],ws);

Au = gs(Au,Qx,Qy);

Au = mask(Au,M);

end
