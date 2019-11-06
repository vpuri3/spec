%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [w] = mass2d(Bmd,Jd,u);

w = ABu(Jd,Jd,u);
w = Bmd .* w;
w = ABu(Jd',Jd',w);

end
