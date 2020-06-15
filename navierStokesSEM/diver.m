%
%	DD
%	(q,dudx+dvdy)
%
function [Du] = diver(ux,uy,Qx2,Qy2,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	[uxdx,~   ] = grad(ux,Drv,Dsv,rxv,ryv,sxv,syv);
	[~   ,uydy] = grad(uy,Drv,Dsv,rxv,ryv,sxv,syv);

	Du = ABu(Jspv',Jrpv',Bv.*(uxdx+uydy));

	Du = gs(Du,Qx2,Qy2);

end
