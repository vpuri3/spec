%
function [egy] = L2(u,B,Qx,Qy);
	vol = sum(sum(B));

	uu = ABu(Qy,Qx,u);
	egy = dot(uu.*uu,B)/vol;
end
