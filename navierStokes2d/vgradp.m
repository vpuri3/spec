% DD'
% (v,-\grad p)
function [px,py] = vgradp(p,Bp,Jrvp,Jsvp,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	JBp = ABu(Jsvp',Jrvp',Bp .* p);

	px = ABu(Isv,Drv',rxv .* JBp) + ABu(Dsv',Irv,sxv .* JBp);
	py = ABu(Isv,Drv',ryv .* JBp) + ABu(Dsv',Irv,syv .* JBp);

end


