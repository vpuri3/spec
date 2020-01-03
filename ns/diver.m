% DD
% (q,dudx+dvdy)
function [Du] = diver(ux,uy,Qx1,Qy1,Qx2,Qy2,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	uux = ABu(Qy1,Qx1,ux);
	uuy = ABu(Qy1,Qx1,uy);

	[uxdx,~   ] = grad(uux,Drv,Dsv,rxv,ryv,sxv,syv);
	[~   ,uydy] = grad(uuy,Drv,Dsv,rxv,ryv,sxv,syv);

	Du = ABu(Jspv',Jrpv',Bv.*(uxdx+uydy));
	Du = ABu(Qy2',Qx2',Du);

end


