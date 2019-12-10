%
%-------------------------------------------------------------------------------
% Preconditioned Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------
function [x] = pcg_pres(b,x0,tol,maxiter...
		   	   ,Bv,Biv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv...
			   ,Rxvx,Ryvx,Rxvy,Ryvy,Sxp,Syp,Lip)
x  = x0;
Ax=pres_E(x,Bv,Biv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv,Rxvx,Ryvx,Rxvy,Ryvy);
ra = b - Ax;
ha  = 0;
hp  = 0;
hpp = 0;
rp  = 0;
rpp = 0;
u   = 0;
k   = 0;

while norm(ra, inf) > tol
	ha = fdm(ra,Sxp,Syp,Lip);
	k = k + 1;
	if (k==maxiter),warning('no conversion.'); return; end;
	hpp = hp;
	rpp = rp;
	hp = ha;
	rp = ra;
	t = dot(rp,hp);
	if k == 1; u = hp;
	else; u = hp + (t / dot(rpp,hpp)) * u;
	end
	Au=pres_E(u,Bv,Biv,Jrpv,Jspv,Irv,Isv,Drv,Dsv,rxv,ryv,sxv,syv,Rxvx,Ryvx,Rxvy,Ryvy);
	a = t / dot(u,Au);
	x = x + a * u;
	ra = rp - a * Au;
end;
['pcg iter:',num2str(k),', residual:',num2str(norm(ra,inf))];

end
