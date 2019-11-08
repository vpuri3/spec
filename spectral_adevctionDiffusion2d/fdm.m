%
function [u] = fdm(b,B,Sr,Ss,D);

u = b ./ B;
u = ABu(Ss',Sr',u);
u = u ./ D;
u = ABu(Ss,Sr,u);

end
