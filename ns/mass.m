%
% input: Bd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [Bu] = mass(u,Bd,M,Qx,Qy,Jr,Js);

uu = mask(u,M);
uu = ABu(Qy,Qx,uu);

if(nargin==5); Jr=[]; Js=[]; end;

ud  = ABu(Js,Jr,uu);
Bud = Bd .* ud;
Bu  = ABu(Js',Jr',Bud);

Bu = ABu(Qy',Qx',Bu);
Bu = mask(Bu,M);

end
