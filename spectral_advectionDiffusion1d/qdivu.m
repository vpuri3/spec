% DD
% (q,dudx+dvdy)
function [q] = qdivu(ux,uy,mskx,msky,Bp,Jrvp,Jsvp,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)
	uux = mask(ux,mskx);
	uuy = mask(uy,msky);

	[uxdx,uxdy] = grad(uux,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);
	[uydx,uydy] = grad(uuy,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv);

	q = Bp .* ABu(Jsvp,Jrvp,uxdx + uydy);
end


