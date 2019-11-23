% pressure project

function [vx,vy,pr] = pres_proj(ux,uy,pr1...
							,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Bpi,Srp,Ssp,Srip,Ssip,Lpi)

	g = -qdivu(ux,uy,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	if(slv==1) % FDM
		pr = fdm(g,Bpi,Srp,Ssp,Srip,Ssip,Lpi);
	end

	[px,py] = vgradp(pr,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	vx = ux + 1/(b0)*RBivx.*px;
	vy = uy + 1/(b0)*RBivy.*py;

	pr = pr + pr1;

end

