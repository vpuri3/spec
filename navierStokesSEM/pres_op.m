%
function [Ep] = pres_op(p,Mvx,Mvy,Qx1,Qy1,Qx2,Qy2...
					   ,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	% DD'
	[px,py] = gradp(p,Qx1,Qy1,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

	% approx hlmhltz inv
	px = mass(px,Biv,Mvx,Qx1,Qy1);
	py = mass(py,Biv,Mvy,Qx1,Qy1);

	% DD
	Ep = diver(px,py,Qx2,Qy2,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);

end
