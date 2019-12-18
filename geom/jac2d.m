%
function [J,Ji,rx,ry,sx,sy] = jac(x,y,Ir,Is,Dr,Ds);

xr = ABu(Is,Dr,x); % dx/dr
xs = ABu(Ds,Ir,x);
yr = ABu(Is,Dr,y);
ys = ABu(Ds,Ir,y);

J  = xr .* ys - xs .* yr;
Ji = 1 ./ J;

rx =  Ji .* ys;
ry = -Ji .* xs;
sx = -Ji .* yr;
sy =  Ji .* xr;

end
