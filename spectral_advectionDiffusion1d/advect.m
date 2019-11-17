%
%   (\vect{v},\vect{c} \cdot \grad(u))
%
function [w] = advect(u,msk,cx,Bmd,J,Dr);

uu = mask(u,msk);

[ux] = grad(uu,Dx);

uxd = J*ux;
cxd = J*cx;

w = Bmd .* (cxd.*uxd);
w = J'*w;

w = mask(w,msk);

end
