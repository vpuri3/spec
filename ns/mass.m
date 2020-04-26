%
% input: B - diag mass mat
%
function [Bu] = mass(u,B,M,Qx,Qy);

Mu = mask(u,M);

if(length(B)==0); Bu=   Mu;
else			  Bu=B.*Mu;
end

Bu = gs(Bu,Qx,Qy);
Bu = mask(Bu,M);

end
