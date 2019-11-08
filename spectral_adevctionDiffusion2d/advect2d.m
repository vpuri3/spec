%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [w] = advect2d(u,cx,cy,Bmd,Ir,Is,Jr,Js,Dr,Ds,rx,ry,sx,sy);

[ux,uy] = grad2d(u,Ir,Is,Dr,Ds,rx,ry,sx,sy);

uxd = ABu(Js,Jr,ux);
uyd = ABu(Js,Jr,uy);
cxd = ABu(Js,Jr,cx);
cyd = ABu(Js,Jr,cy);

w = cxd.*uxd + cyd.*uyd;
w = Bmd.*w;
w = ABu(Js',Jr',w);

end
