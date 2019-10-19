%
%     -\del^2 u
%
function [v] = laplace(Dm1,Jd,Bmd,rxmd,rymd,sxmd,symd,u);

wr = Bmd.*ABu(Jd,Jd*Dm1,u);
ws = Bmd.*ABu(Jd*Dm1,Jd,u);

G11 = rxmd.*rxmd + rymd.*rymd;
G12 = rxmd.*sxmd + rymd.*symd;
G22 = sxmd.*sxmd + symd.*symd;

vr = G11.*wr + G12.*ws;
vs = G12.*wr + G22.*ws;

v = ABu(Jd',Dm1*Jd',vr) + ABu(Dm1*Jd',Jd',vs);

end
