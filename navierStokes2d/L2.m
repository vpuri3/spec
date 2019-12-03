%
function [egy] = L2(u,B);
	vol = sum(sum(B));
	egy = dot(u.*u,B)/vol;
end
