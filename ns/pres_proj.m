%
% pressure projection operator
%
function [vx,vy,pr] = pres_proj(ux,uy,pr,b0,Mvx,Mvy,Qx1,Qy1,Qx2,Qy2...
							,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	g = -diver(ux,uy,Qx2,Qy2,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	delp = b0 * pcg_pres(g,g,1e-8,1e4...
	   	   ,Mvx,Mvy,Qx1,Qy1,Qx2,Qy2,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	[px,py] = vgradp(delp,Qx1,Qy1,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	px = mass(px,Biv/b0,Mvx,Qx1,Qy1);
	py = mass(py,Biv/b0,Mvy,Qx1,Qy1);

	vx = ux + px;
	vy = uy + py;

	pr = pr + delp;

end
