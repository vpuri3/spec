%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [w] = mass2d(u,msk,Bmd,Jr,Js);

w = mask(u,msk);

w = ABu(Js,Jr,u);
w = Bmd .* w;
w = ABu(Js',Jr',w);

w = mask(w,msk);

end
