%
% viscous solve
%
function [u] = visc_slv(b,Sr,Ss,Li,slv...
					   ,Rx,Ry,visc,b0,B,Ir,Is,Dr,Ds,g11,g12,g22)
	if(slv==0) % CG

		u = pcg_visc(b,0*b,1e-8,1e3...
				   ,Rx,Ry,visc,b0,B,Ir,Is,Dr,Ds,g11,g12,g22...
				   ,Sr,Ss,Li);

	elseif(slv==1) % FDM

		u = fdm(b,Sr,Ss,Li);

	end	
end
