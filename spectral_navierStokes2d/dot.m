%
function [a] = dot(v,u);
	a = v .* u;
	a = sum(sum(a));
end
