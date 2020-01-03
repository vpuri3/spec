%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [Cu] = advect(u,cx,cy,M,Qx,Qy,Bmd,Jr,Js,Dr,Ds,rx,ry,sx,sy);

uu  = mask(u,M);

uu  = ABu(Qy,Qx,uu);
cxl = ABu(Qy,Qx,cx);
cyl = ABu(Qy,Qx,cy);

[ux,uy] = grad(uu,Dr,Ds,rx,ry,sx,sy);

uxd = ABu(Js,Jr,ux);
uyd = ABu(Js,Jr,uy);
cxd = ABu(Js,Jr,cxl);
cyd = ABu(Js,Jr,cyl);

Cud = cxd.*uxd + cyd.*uyd;
Cud = Bmd.*Cud;
Cu  = ABu(Js',Jr',Cud);

Cu = ABu(Qy',Qx',Cu);
Cu = mask(Cu,M);

end
