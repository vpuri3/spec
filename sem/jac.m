%
function [J,Jinv,rx,ry,sx,sy] = jac(x,y,I,D);

xr = ABu(D,I,x); % dx/dr
xs = ABu(I,D,x);
yr = ABu(D,I,y);
ys = ABu(I,D,y);

J    = xr.*ys - xs*yr;
Jinv = 1 ./ J;

rx = Jinv.*ys;
ry =-Jinv.*xs;
sx =-Jinv.*yr;
sy = Jinv.*xr;

end
