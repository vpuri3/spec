% DD'
% (v,-\grad p)
function [px,py] = vgradp(p,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	BJp = Bv .* ABu(Jspv,Jrpv,p);

	px = ABu([],Drv',rxv .* BJp) + ABu(Dsv',[],sxv .* BJp);
	py = ABu([],Drv',ryv .* BJp) + ABu(Dsv',[],syv .* BJp);

end


