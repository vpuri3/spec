%
%     (v,-\del^2 u)
%
function [w] = laplace2d(u,Dm1,Jd,G11,G12,G22);

urd = ABu(Jd,Jd*Dm1,u);
usd = ABu(Jd*Dm1,Jd,u);

wr = G11.*urd + G12.*usd;
ws = G12.*urd + G22.*usd;

w = ABu(Jd',Dm1'*Jd',wr) + ABu(Dm1'*Jd',Jd',ws);

end

