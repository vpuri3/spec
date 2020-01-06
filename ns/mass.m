%
% input: Bd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [Bu] = mass(u,Bd,M,Qx,Qy,Jr,Js);

uu = mask(u,M);

if(nargin==5); Jr=[]; Js=[]; end;

ud  = ABu(Js,Jr,uu);
Bud = Bd .* ud;
Bu  = ABu(Js',Jr',Bud);

Bu = ABu(Qy',Qx',Bu); % gather
Bu = ABu(Qy ,Qx ,Bu); % scatter

Bu = mask(Bu,M);

end
