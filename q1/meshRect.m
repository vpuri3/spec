% rectangle mesh
ne = 4; EE = zeros(ne,4);
% EE(ei,:) is the list of vertices associated with element ei.
EE(1,:) = [1,2,5,4];
EE(2,:) = [2,3,6,5];
EE(3,:) = [4,5,8,7];
EE(4,:) = [5,6,9,8];
YY = [0,0,0,0.5,0.5,0.5,1,1,1]';
XX = [0,0.5,1,0,0.5,1,0,0.5,1]';
lx = 1/sqrt(ne); ly=lx;