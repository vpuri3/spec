%
%	DD'
%	(v,-\grad p)
%
function [px,py] = vgradp(p,Qx1,Qy1,Bv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv)

	pp = Bv .* ABu(Jspv,Jrpv,p);

	px = ABu([],Drv',rxv .* pp) + ABu(Dsv',[],sxv .* pp);
	py = ABu([],Drv',ryv .* pp) + ABu(Dsv',[],syv .* pp);

	px = ABu(Qy1',Qx1',px); % gather
	py = ABu(Qy1',Qx1',py);

	px = ABu(Qy1,Qx1,px); % scatter
	py = ABu(Qy1,Qx1,py);

end

