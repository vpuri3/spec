%
% viscous solve
%
function [u] = visc_slv(b,visc,b0,M,Qx,Qy,Bi,B,Dr,Ds,g11,g12,g22)

	u = pcg_visc(b,0*b,1e-8,1e3...
			   ,M,Qx,Qy,visc,b0,B,Bi,Dr,Ds,g11,g12,g22);
end
