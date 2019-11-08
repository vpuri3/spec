%
function [u] = fdm(b,B,Sr,Ss,Rx,Ry,D);

u = b ./ B;
u = ABu(Ry ,Rx ,u);
u = ABu(Ss',Sr',u);
u = u ./ D;
u = ABu(Ss ,Sr ,u);
u = ABu(Ry',Rx',u);

end
