%
function [w] = extend(u,Rx,Ry);

w = Rx' * u * Ry;
