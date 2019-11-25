%
% pressure project
%
function [vx,vy,pr] = pres_proj(ux,uy,pr1...
							,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy,slv...
							,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Bip,Srp,Ssp,Srip,Ssip,Lip)

	g = -qdivu(ux,uy,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	if(slv==1) % FDM
		delp = (1/b0)*fdm(g,Bip,Srp,Ssp,Srip,Ssip,Lip);
	end

	[px,py] = vgradp(delp,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	RBpx = ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,(1/b0) * Biv .* px);
	RBpy = ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,(1/b0) * Biv .* py);

	vx = ux + RBpx;
	vy = uy + RBpy;

	pr = delp + pr1;

end

