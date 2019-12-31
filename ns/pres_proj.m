%
% pressure projection operator
%
function [vx,vy,pr] = pres_proj(ux,uy,pr,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy...
							,Bp,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Bv)

	g = -diver(ux,uy,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	delp = b0 * pcg_pres(g,0*g,1e-8,1e3...
	   	   ,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv...
		   ,Rxvx,Ryvx,Rxvy,Ryvy);

	[px,py] = vgradp(delp,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	dpvx = (1/b0) * Biv .* ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	dpvy = (1/b0) * Biv .* ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	dpvx = ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,dpvx);
	dpvy = ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,dpvy);

	vx = ux + dpvx;
	vy = uy + dpvy;

	pr = pr + delp;

end
