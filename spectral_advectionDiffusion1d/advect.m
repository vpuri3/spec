%
%   (v,c*u_x)
%
function [Cu] = advect(u,c,Bd,D,J);

uxd = J*D*u;
cd  = J*c;
Cud = Bd .* (cd.*uxd);

Cu = J'*Cud;

end
