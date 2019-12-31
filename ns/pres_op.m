%
function [Ep] = pres_op(p,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv...
					  ,Rxvx,Ryvx,Rxvy,Ryvy)

	% DD'
	[px,py] = vgradp(p,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	% QQ
	px = ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	py = ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	px = Biv .* px;	
	py = Biv .* py;

	px = ABu(Ryvx'*Ryvx,Rxvx'*Rxvx,px);
	py = ABu(Ryvy'*Ryvy,Rxvy'*Rxvy,py);

	% DD
	Ep = diver(px,py,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

end

