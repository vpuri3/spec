%
%	DD'
%	(v,-\grad p)
%
function [px,py] = vgradp(p,Qx1,Qy1,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	pp = Bv .* ABu(Jspv,Jrpv,p);

	px = ABu([],Drv',rxv .* pp) + ABu(Dsv',[],sxv .* pp);
	py = ABu([],Drv',ryv .* pp) + ABu(Dsv',[],syv .* pp);

	px = gs(px,Qx1,Qy1);
	py = gs(py,Qx1,Qy1);

end

