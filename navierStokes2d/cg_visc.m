%
%-------------------------------------------------------------------------------
% Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------
function [x] = cg_visc(b,x0,tol,maxiter...
				      ,Rx,Ry,visc,b0,B,Ir,Is,Dr,Ds,g11,g12,g22...
				   	  ,Sr,Ss,Li)
	x  = x0;
	Ax = ABu(Ry,Rx,hlmhltz(ABu(Ry',Rx',x),visc,b0,B,Ir,Is,Dr,Ds,g11,g12,g22));
	r  = b - Ax;
	rsqold=dot(r,r);
	
	if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;
	
	p=r;
	for k=1:maxiter
		Ap = ABu(Ry,Rx,hlmhltz(ABu(Ry',Rx',p),visc,b0,B,Ir,Is,Dr,Ds,g11,g12,g22));
		al = rsqold / dot(p,Ap);
		x  = x + al*p;
		r  = r - al*Ap;
		rsqnew=dot(r,r); if(sqrt(rsqnew) < tol); break; end;
		be = rsqnew / rsqold;
		p  = r + be*p;
		rsqold = rsqnew;
	end
	['cg iter:',num2str(k),', residual:',num2str(sqrt(rsqnew))]
end

