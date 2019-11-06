%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [w] = advec(u,cx,cy,Dm1,J1d,Bmd,rxmd,rymd,sxmd,symd);

[ux,uy] = grad2d(u,Ir,Is,Dr,Ds,rxm1,rym1,sxm1,sym1);

uxd = ABu(Js,Jr,ux);
uyd = ABu(Js,Jr,uy);
cxd = ABu(Js,Jr,cx);
cyd = ABu(Js,Jr,cy);

w = cxd.*uxd + cyd.*uyd;
w = Bmd.*w;
w = ABu(Js',Jr',w);

end
