%
%     (v,-\del^2 u)
%
function [Au] = lapl(u,B,D);

Au = D'*B*D*u;

end

