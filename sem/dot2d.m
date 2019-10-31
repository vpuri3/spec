%
function [a] = dot2d(v,u,Bdiag);
	a = v .* Bdiag .* u;
	a = sum(sum(a));
end
