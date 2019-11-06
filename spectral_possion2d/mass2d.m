%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [w] = mass2d(Bmd,Jr,Js,u);

w = ABu(Js,Jr,u);
w = Bmd .* w;
w = ABu(Js',Jr',w);

end
