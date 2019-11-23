%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [Cu] = advect(u,cx,cy,Bmd,Ir,Is,Jr,Js,Dr,Ds,rx,ry,sx,sy);

[ux,uy] = grad(u,Ir,Is,Dr,Ds,rx,ry,sx,sy);

uxd = ABu(Js,Jr,ux);
uyd = ABu(Js,Jr,uy);
cxd = ABu(Js,Jr,cx);
cyd = ABu(Js,Jr,cy);

Cud = cxd.*uxd + cyd.*uyd;
Cud = Bmd.*Cud;
Cu  = ABu(Js',Jr',Cud);

end
