%
% (q,dudx+dvdy)
function [DDu] = diver(ux,uy,Bp,Jrvp,Jsvp,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	[uxdx,uxdy] = grad(ux,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);
	[uydx,uydy] = grad(uy,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	DDu = Bp .* ABu(Jsvp,Jrvp,uxdx + uydy);
end


