%
% (q,dudx+dvdy)
function [DDu] = diver(ux,uy,Bv,Jrvp,Jsvp,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	[uxdx,uxdy] = grad(ux,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);
	[uydx,uydy] = grad(uy,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	DDu = ABu(Jsvp,Jrvp,Bv.*(uxdx+uydy));
end


