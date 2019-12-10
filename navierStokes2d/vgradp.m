% DD'
% (v,-\grad p)
function [px,py] = vgradp(p,Bv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	BJp = Bv .* ABu(Jspv,Jrpv,p);

	px = ABu(Isv,Drv',rxv .* BJp) + ABu(Dsv',Irv,sxv .* BJp);
	py = ABu(Isv,Drv',ryv .* BJp) + ABu(Dsv',Irv,syv .* BJp);

end


