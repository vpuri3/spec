%
function [J,Jinv,rx,ry,sx,sy] = jac(x,y,Ir,Is,Dr,Ds);

xr = ABu(Is,Dr,x); % dx/dr
xs = ABu(Ds,Ir,x);
yr = ABu(Is,Dr,y);
ys = ABu(Ds,Ir,y);

J    = xr.*ys - xs.*yr;
Jinv = 1 ./ J;

rx =  Jinv.*ys;
ry = -Jinv.*xs;
sx = -Jinv.*yr;
sy =  Jinv.*xr;

end
