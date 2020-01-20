%
%	DD'
%	(v,-\grad p) = (\grad v,p)
%
function [px,py] = gradp(p,Qx1,Qy1,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	Jp = ABu(Jspv,Jrpv,p);
	Bp = mass(Jp,Bv,[],[],[]);

	px = ABu([],Drv',rxv .* Bp) + ABu(Dsv',[],sxv .* Bp);
	py = ABu([],Drv',ryv .* Bp) + ABu(Dsv',[],syv .* Bp);

	px = gs(px,Qx1,Qy1);
	py = gs(py,Qx1,Qy1);

end

