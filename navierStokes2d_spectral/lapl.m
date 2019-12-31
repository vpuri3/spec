%
%     (v,-\del^2 u)
%
function [Au] = lapl(u,Ir,Is,Dr,Ds,G11,G12,G22)

ur = ABu(Is,Ir*Dr,u);
us = ABu(Is*Ds,Ir,u);

wr = G11.*ur + G12.*us;
ws = G12.*ur + G22.*us;

Au = ABu(Is',Dr'*Ir',wr) + ABu(Ds'*Is',Ir',ws);

end

