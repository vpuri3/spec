%
function [x] = cg2d(A,b,x0,tol,maxiter);

x =x0;
r0=b-Afun(x);
mag=sum(sum(Bm1 .* r0)); if(mag < tol); x=x0; return; end;

maxiter = min(maxiter,length(x0);
p0=r0;
for k=1:maxiter
	
	al = dot2d(p1,b   ,Bm1) / dot2d(p1,A*pm1,Bm1); % alpha
	be = dot2d(p1,A*b),Bm1) / dot2d(p1,A*pm1,Bm1); % beta

	p = p - be*r;
	x = x + al*p;
	r = r - al*p;

end

end
