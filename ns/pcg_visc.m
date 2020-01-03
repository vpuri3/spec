%
%-------------------------------------------------------------------------------
% Preconditioned Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------
function [x] = pcg_visc(b,visc,b0,M,Qx,Qy...
					   ,B,Bi,Dr,Ds,g11,g12,g22,tol,maxiter);
x  = b;
Ax = hlmhltz(x,visc,b0,M,Qx,Qy,B,Dr,Ds,g11,g12,g22);
ra = b - Ax;
ha  = 0;
hp  = 0;
hpp = 0;
rp  = 0;
rpp = 0;
u   = 0;
k   = 0;

while norm(ra, inf) > tol
	ha = mass(ra,Bi/b0,M,Qx,Qy);
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
	Au = hlmhltz(u,visc,b0,M,Qx,Qy,B,Dr,Ds,g11,g12,g22);
	a = t / dot(u,Au);
	x = x + a * u;
	ra = rp - a * Au;
end;
['pcg visc:',num2str(k),', residual:',num2str(norm(ra,inf))];

end
