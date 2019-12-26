%
% viscous solve
%
function [u] = visc_slv(b...
					,Bi,Rx,Ry,visc,b0,B,Ir,Is,Dr,Ds,g11,g12,g22)

	u = pcg_visc(b,0*b,1e-8,1e3...
			   ,Rx,Ry,visc,b0,B,Bi,Ir,Is,Dr,Ds,g11,g12,g22);
end
