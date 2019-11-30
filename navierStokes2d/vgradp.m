% DD'
% (v,-\grad p)
function [px,py] = vgradp(p,Bv,Jrvp,Jsvp,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	JBp = Bv .* ABu(Jsvp',Jrvp',p);

	px = ABu(Isv,Drv',rxv .* JBp) + ABu(Dsv',Irv,sxv .* JBp);
	py = ABu(Isv,Drv',ryv .* JBp) + ABu(Dsv',Irv,syv .* JBp);

end


