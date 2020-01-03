% DD'
% (v,-\grad p)
function [px,py] = vgradp(p,Qx1,Qy1,Qx2,Qy2,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	pp = ABu(Qy2,Qx2,p);

	Bp = Bv .* ABu(Jspv,Jrpv,pp);

	px = ABu([],Drv',rxv .* Bp) + ABu(Dsv',[],sxv .* Bp);
	py = ABu([],Drv',ryv .* Bp) + ABu(Dsv',[],syv .* Bp);

	px = ABu(Qy1',Qx1',px);
	py = ABu(Qy1',Qx1',py);

end


