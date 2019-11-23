%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [Bu] = mass(u,Bmd,Jr,Js);

ud  = ABu(Js,Jr,u);
Bud = Bmd .* ud;
Bu  = ABu(Js',Jr',Bud);

end
