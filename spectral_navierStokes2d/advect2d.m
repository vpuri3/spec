%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [w] = advect2d(u,msk,cx,cy,Bmd,Ir,Is,Jr,Js,Dr,Ds,rx,ry,sx,sy);

uu = mask(u,msk);

[ux,uy] = grad2d(uu,Ir,Is,Dr,Ds,rx,ry,sx,sy);

uxd = ABu(Js,Jr,ux);
uyd = ABu(Js,Jr,uy);
cxd = ABu(Js,Jr,cx);
cyd = ABu(Js,Jr,cy);

w = cxd.*uxd + cyd.*uyd;
w = Bmd.*w;
w = ABu(Js',Jr',w);

w = msk(w,msk);

end
