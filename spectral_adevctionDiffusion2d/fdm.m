%
function [u] = fdm(b,B,Sr,Ss,Sri,Ssi,Rx,Ry,D);

u = b ./ B;
u = ABu(Ry ,Rx ,u);
u = ABu(Ssi,Sri,u);
u = u ./ D;
u = ABu(Ss ,Sr ,u);
u = ABu(Ry',Rx',u);

end
