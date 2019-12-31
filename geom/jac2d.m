%
function [J,Ji,rx,ry,sx,sy] = jac(x,y,Dr,Ds);

xr = ABu([],Dr,x); % dx/dr
xs = ABu(Ds,[],x);
yr = ABu([],Dr,y);
ys = ABu(Ds,[],y);

J  = xr .* ys - xs .* yr;
Ji = 1 ./ J;

rx =  Ji .* ys;
ry = -Ji .* xs;
sx = -Ji .* yr;
sy =  Ji .* xr;

end
