%
% pressure projection operator
%
function [vx,vy,pr] = pres_proj(ux,uy,pr,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy,slv...
							,Bp,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Sxp,Syp,Lip,Bv)

	g = -diver(ux,uy,Bv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	if(slv==1) % FDM
		delp = b0 * fdm(g,Sxp,Syp,Lip);
	else
		delp = b0 * pcg_pres(g,0*g,1e-8,1e3...
		   	   ,Bv,Biv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
			   ,Rxvx,Ryvx,Rxvy,Ryvy,Sxp,Syp,Lip);
	end

	[px,py] = vgradp(delp,Bv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	dpvx = (1/b0) * Biv .* ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	dpvy = (1/b0) * Biv .* ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	dpvx = ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,dpvx);
	dpvy = ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,dpvy);

	vx = ux + dpvx;
	vy = uy + dpvy;

	pr = pr + delp;

end
