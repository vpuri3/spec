%
function [J,Jinv,rx,ry,sx,sy] = jac(x,y,I,D);

xr = ABu(I,D,x); % dx/dr
xs = ABu(D,I,x);
yr = ABu(I,D,y);
ys = ABu(D,I,y);

J    = xr.*ys - xs*yr;
Jinv = 1 ./ J;

rx = Jinv.*ys;
ry =-Jinv.*xs;
sx =-Jinv.*yr;
sy = Jinv.*xr;

end
