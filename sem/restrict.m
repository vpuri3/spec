%
function [w] = restrict(u,Rx,Ry);

w = Rx * u * Ry';
