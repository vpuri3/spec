%
% pressure project
%
function [vx,vy,pr] = pres_proj(ux,uy,pr1...
							,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy,slv...
							,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Srp,Ssp,Lip)

	g = -diver(ux,uy,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	if(slv==1) % FDM
		delp = fdm(g,Srp,Ssp,Lip);
	end

	[px,py] = vgradp(delp,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	px = (1/b0) * Biv .* ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	py = (1/b0) * Biv .* ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	vx = ux + px;
	vy = uy + py;

	pr = pr1 + delp;

end
