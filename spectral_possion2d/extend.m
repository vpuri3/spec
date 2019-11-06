%
function [w] = extend(u,Rx,Ry);

w = ABu(Ry',Rx',u);
