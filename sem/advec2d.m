%
%   \vect{c} \cdot \grad(u)
%
function [v] = advec(cx,cy,Dm1,Jd,Bmd,rxmd,rymd,sxmd,symd,u,msk);

[ux,uy] = grad2d(Jd,Jd*Dm1,rxmd,rymd,sxmd,symd,u);

v = cx.*ux + cy.*uy;

v = Bmd.*v;

v = ABu(Jd',Jd',u);

v = v .* msk;

end
