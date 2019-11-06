%
% from
% https://en.wikipedia.org/wiki/Conjugate_gradient_method
%
function [x,k,rsqnew] = cg_possion2d(b,x0,tol,maxiter,Jr,Js,Dr,Ds,G11,G12,G22,msk);
x = x0;
r = b - msk.*laplace2d(x,Jr,Js,Dr,Ds,G11,G12,G22); % r = b - A*x;
rsqold=dot2d(r,r);

if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;

p=r;
for k=1:maxiter
	Ap = msk.*laplace2d(p,Jr,Js,Dr,Ds,G11,G12,G22); % Ap = A*p;
	al = rsqold / dot2d(p,Ap);
	x  = x + al*p;
	r  = r - al*Ap;
	rsqnew=dot2d(r,r); if(sqrt(rsqnew) < tol); return; end;
	be = rsqnew / rsqold;
	p  = r + be*p;
	rsqold = rsqnew;
end

end
