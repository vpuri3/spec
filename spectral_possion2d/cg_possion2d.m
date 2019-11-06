%
% from
% https://en.wikipedia.org/wiki/Conjugate_gradient_method
%
function [x,k,rsqnew] = cg_possion2d(b,x0,tol,maxiter,Dm1,Jd,G11,G12,G22,msk);
x = x0;
r = b - laplace2d(x,Dm1,Jd,G11,G12,G22,msk); % r = b - A*x;
rsqold=dot2d(r,r);

if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;

p=r;
for k=1:maxiter
	Ap = laplace2d(p,Dm1,Jd,G11,G12,G22,msk); % Ap = A*p;
	al = rsqold / dot2d(p,Ap);
	x  = x + al*p;
	r  = r - al*Ap;
	rsqnew=dot2d(r,r); if(sqrt(rsqnew) < tol); return; end;
	be = rsqnew / rsqold;
	p  = r + be*p;
	rsqold = rsqnew;
end

end
