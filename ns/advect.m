%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [Cu] = advect(u,cx,cy,M,Qx,Qy,Bmd,Jr,Js,Dr,Ds,rx,ry,sx,sy);

uu  = mask(u,M);

[ux,uy] = grad(uu,Dr,Ds,rx,ry,sx,sy);

uxd = ABu(Js,Jr,ux);
uyd = ABu(Js,Jr,uy);
cxd = ABu(Js,Jr,cx);
cyd = ABu(Js,Jr,cy);

Cud = cxd.*uxd + cyd.*uyd;
Cud = Bmd.*Cud;
Cu  = ABu(Js',Jr',Cud);

Cu = gs(Cu,Qx,Qy);
Cu = mask(Cu,M);

end
