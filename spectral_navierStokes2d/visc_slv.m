%
% viscous solve
%
function [u] = visc_slv(b,Sr,Ss,Li,slv)

	if(slv==0) % CG
	
		%u = cg_visc(b,msk,0*b,1e-8,1e3);
	
	elseif(slv==1) % FDM
	
		u = fdm(b,Sr,Ss,Li);
	end	
end

%-------------------------------------------------------------------------------
% Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------
%function [x,k,rsqnew] = cg_visc(b,msk,x0,tol,maxiter);
%	x = x0;
%	r = b - hmhltz(x,msk); % r = b - Ax
%	rsqold=dot(r,r);
%	
%	if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;
%	
%	p=r;
%	for k=1:maxiter
%		Ap = hmhltz(p,msk); % Ap = A*p
%		al = rsqold / dot(p,Ap);
%		x  = x + al*p;
%		r  = r - al*Ap;
%		rsqnew=dot(r,r); if(sqrt(rsqnew) < tol); return; end;
%		be = rsqnew / rsqold;
%		p  = r + be*p;
%		rsqold = rsqnew;
%	end
%	['cg iter:',num2str(k),', residual:',num2str(sqrt(rsqnew))]
%end


