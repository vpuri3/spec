%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [Bu] = mass(u,Bmd,M,Qx,Qy,Jr,Js);

ul = ABu(Qy,Qx,u);
ul = mask(ul,M);

if(nargin==5); Bu = Bmd .* ul; return; end;

ud  = ABu(Js,Jr,ul);
Bud = Bmd .* ud;
Bu  = ABu(Js',Jr',Bud);

Bu = mask(Bu,M);
Bu = ABu(Qy',Qx',Bu);

end
