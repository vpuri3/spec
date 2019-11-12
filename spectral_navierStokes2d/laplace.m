%
%     (v,-\del^2 u)
%
function [w] = laplace(u,msk,Jr,Js,Dr,Ds,G11,G12,G22);

uu = mask(u,msk);

urd = ABu(Js,Jr*Dr,uu);
usd = ABu(Js*Ds,Jr,uu);

wr = G11.*urd + G12.*usd;
ws = G12.*urd + G22.*usd;

w = ABu(Js',Dr'*Jr',wr) + ABu(Ds'*Js',Jr',ws);

w = mask(w,msk);

end

