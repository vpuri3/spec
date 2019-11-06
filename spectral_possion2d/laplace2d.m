%
%     (v,-\del^2 u)
%
function [w] = laplace2d(u,Dm1,Jd,Bmd,rxmd,rymd,sxmd,symd,msk);

ur = Bmd.*ABu(Jd,Jd*Dm1,u); % error in moving bw q-nodes and v-nodes
us = Bmd.*ABu(Jd*Dm1,Jd,u);

G11 = rxmd.*rxmd + rymd.*rymd;
G12 = rxmd.*sxmd + rymd.*symd;
G22 = sxmd.*sxmd + symd.*symd;

wr = G11.*ur + G12.*us;
ws = G12.*ur + G22.*us;

w = ABu(Jd',Dm1'*Jd',wr) + ABu(Dm1'*Jd',Jd',ws); % back to velocity nodes

w = msk .* w;

end
