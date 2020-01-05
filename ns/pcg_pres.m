%
%-------------------------------------------------------------------------------
% Preconditioned Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------
function [x] = pcg_pres(b,x0,tol,maxiter...
		   	   ,Mvx,Mvy,Qx1,Qy1,Qx2,Qy2,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);
x  = x0;
Ax=pres_op(x,Mvx,Mvy,Qx1,Qy1,Qx2,Qy2,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);
ra = b - Ax;
ha  = 0;
hp  = 0;
hpp = 0;
rp  = 0;
rpp = 0;
u   = 0;
k   = 0;

while norm(ra,inf) > tol
	ha = ra;
	k = k + 1;
	if (k==maxiter),warning('no conversion.'); return; end;
	hpp = hp;
	rpp = rp;
	hp = ha;
	rp = ra;
	t = dot(rp,hp);
	if k == 1; u = hp; else; u = hp + (t / dot(rpp,hpp)) * u; end
	Au=pres_op(u,Mvx,Mvy,Qx1,Qy1,Qx2,Qy2,Bv,Biv,Jrpv,Jspv,Drv,Dsv,rxv,ryv,sxv,syv);
	a = t / dot(u,Au);
	x = x + a * u;
	ra = rp - a * Au;
end;
['pcg pres iter: ',num2str(k),', residual: ',num2str(norm(ra,inf))];

end
