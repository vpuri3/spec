%
%     (v,-\del^2 u)
%
function [Au] = lapl(u,Jr,Js,Dr,Ds,G11,G12,G22);

urd = ABu(Js,Jr*Dr,u);
usd = ABu(Js*Ds,Jr,u);

wr = G11.*urd + G12.*usd;
ws = G12.*urd + G22.*usd;

Au = ABu(Js',Dr'*Jr',wr) + ABu(Ds'*Js',Jr',ws);

end

