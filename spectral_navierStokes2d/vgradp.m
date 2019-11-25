% DD'
% (v,-\grad p)
function [qx,qy] = vgradp(q,Bp,Jrvp,Jsvp,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	JBq = ABu(Jsvp',Jrvp',Bp .* q);

	qx = ABu(Isv,Drv',rxv .* JBq) + ABu(Dsv',Irv,sxv .* JBq);
	qy = ABu(Isv,Drv',ryv .* JBq) + ABu(Dsv',Irv,syv .* JBq);

end


