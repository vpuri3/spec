% DD
% (q,dudx+dvdy)
function [DDu] = diver(ux,uy,Bv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	[uxdx,uxdy] = grad(ux,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);
	[uydx,uydy] = grad(uy,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	DDu = ABu(Jspv',Jrpv',Bv.*(uxdx+uydy));
end


