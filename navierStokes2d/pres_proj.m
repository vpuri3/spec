%
% pressure project
%
function [vx,vy,pr] = pres_proj(ux,uy,pr1,b0,Biv,Rxvx,Ryvx,Rxvy,Ryvy,slv...
							,Bp,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
							,Sxp,Syp,Lip,Bv)

	g = -diver(ux,uy,Bv,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);
	['in pres proj']
	mesh(g),title('g'),pause

	if(slv==1) % FDM
		delp = b0 * fdm(g,Sxp,Syp,Lip);
	end

	%delp = delp - dot(Bp,delp)/dot(1+0*Bp,Bp);
	mesh(delp),title('delp'),pause

	[px,py] = vgradp(delp,Bv,Jr,Js,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);
	mesh(px),title('px'),pause
	mesh(py),title('py'),pause

	dpvx = (1/b0) * Biv .* ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	dpvy = (1/b0) * Biv .* ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);
	mesh(dpvx),title('dpvx'),pause
	mesh(dpvy),title('dpvy'),pause

	vx = ux + dpvx;
	vy = uy + dpvy;

	pr = pr1 + delp;

end
