%
% pressure project
%
function [vx,vy,pr] = pres_proj(ux,uy,pr1,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy,slv...
							,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Sxp,Syp,Lip,Bv)

	g = -diver(ux,uy,Bv,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	if(slv==1) % FDM
		delp = b0 * fdm(g,Sxp,Syp,Lip);
	end
	
	delp = delp - dot(Bp,delp)/dot(1+0*Bp,Bp);

	[px,py] = vgradp(delp,Bv,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	px = (1/b0) * Biv .* ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	py = (1/b0) * Biv .* ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	vx = ux + px;
	vy = uy + py;

	pr = pr1 + delp;

end
