%
function [a] = dot2d(v,u);
	a = v .* u;
	a = sum(sum(a));
end
