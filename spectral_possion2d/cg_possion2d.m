%
function [x] = cg_possion2d(A,b,x0,tol,maxiter);

x=x0;
r=b-A*x;
rmg=sqrt(dot2d(r,r)); if(mag < tol); end; end;

maxiter = min(maxiter,length(x0);
p=r;
for k=1:maxiter
	al = rmg*rmg / dot2d(p,A*p);
	x  = x + al*p;
	r  = b - A *x;
	be = dot2d(r,r) / (rmg*rmg);
	rmg=sqrt(dot2d(r,r)); if(mag < tol); break; end;
	p  = r + be*p;
end

end
