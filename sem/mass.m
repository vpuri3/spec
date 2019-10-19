%
% input: Bmd - diag mass mat on quadrature nodes
%        Jd  - interpolation matrix from evaluation
%              to quadrature nodes
function [v] = mass(Bmd,Jd,u);

v = ABu(Jd,Jd,u);
v = Bmd.*v;
v = ABu(Jd',Jd',v);

end
