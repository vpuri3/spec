% DD
% (q,dudx+dvdy)
function [DDu] = diver(ux,uy,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	[uxdx,uxdy] = grad(ux,Drv,Dsv,rxv,ryv,sxv,syv);
	[uydx,uydy] = grad(uy,Drv,Dsv,rxv,ryv,sxv,syv);

	DDu = ABu(Jspv',Jrpv',Bv.*(uxdx+uydy));
end


