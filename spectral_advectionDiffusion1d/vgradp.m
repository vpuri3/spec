% DD'
% (v,-\grad p) = [ (vxdx,p), (vydy,p) ]
function [qx,qy] = vgradp(q,mskx,msky,Bp,Jrvp,Jsvp...
						 ,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv)

	JBq = ABu(Jsvp',Jrvp',Bp .* q);

	qx = ABu(Isv,Drv',rxv.*JBq) + ABu(Dsv',Irv,sxv.*JBq);
	qy = ABu(Isv,Drv',ryv.*JBq) + ABu(Dsv',Irv,syv.*JBq);

	qx = mask(qx,mskx);
	qy = mask(qy,msky);

end


