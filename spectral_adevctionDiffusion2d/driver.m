%======================================================================
%
%	Driver function for advection diffusion equation
%
%	du/dt + \vect{c}\dot\grad{u}  = f + nu*\del^2 u + q*u 
%
%   + Dirichlet/Neumann BC
%
%======================================================================
function driver
%

clf; format compact; format shorte;

nx1 = 16;
ny1 = 16;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);

[zrm1,wrm1] = zwgll(nx1-1);
[zsm1,wsm1] = zwgll(ny1-1);
[zrmd,wrmd] = zwgll(nxd-1);
[zsmd,wsmd] = zwgll(nyd-1);

Drm1 = dhat(zrm1);
Dsm1 = dhat(zsm1);
Drmd = dhat(zrmd);
Dsmd = dhat(zsmd);

Irm1 = eye(nx1);
Ism1 = eye(ny1);
Irmd = eye(nxd);
Ismd = eye(nyd);

Jr1d = interp_mat(zrmd,zrm1); % nx1 to nxd
Js1d = interp_mat(zsmd,zsm1);

%----------------------------------------------------------------------
% deform geometry

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrm1,zsm1);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrmd,zsmd);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xmd,ymd] = ndgrid(zrmd,zsmd);

%----------------------------------------------------------------------
% data

slv=1;                           % solver --> 0: CG, 1: FDM, 2: direct

nu= 1e-1;                        % viscosity
u = sin(pi*xm1).*sin(pi*ym1);    % initial condition
cx= 0+0*xm1;                     % advecting field
cy= 0+0*xm1;
f = sin(pi*xm1).*sin(pi*ym1);    % forcing
f = 1+0*xm1;

ub = u; 					     % dirichlet BC

%Rx = Irm1(1:end-1,:); Rx(end,1)=1;  % periodic
Rx = Irm1(2:end-1,:);                % dir-dir
Ry = Ism1(2:end-1,:);                % dir-dir
% time
T   = 0; % steady state solve if T == 0
CFL = 0.5;

%----------------------------------------------------------------------
% setup

% time stepper
dx = min(min(xm1));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

% BDFk-EXTk setup
a = zeros(3,1);
b = zeros(4,1);

% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% mass
Bm1 = Jm1.*(wrm1*wsm1');
Bmd = Jmd.*(wrmd*wsmd');

% mask
msk = diag(Rx'*Rx) * diag(Ry'*Ry)';

% laplace operator setup
G11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
G12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
G22 = Bmd .* (sxmd.*sxmd + symd.*symd);

if(slv==1) % fast diagonalization setup

	Lx = max(max(xm1))-min(min(xm1));
	Ly = max(max(ym1))-min(min(ym1));
	
	Br = (Lx/2)*diag(wrm1);
	Bs = (Ly/2)*diag(wsm1);
	Dr = (2/Lx)*Drm1;
	Ds = (2/Ly)*Dsm1;
	Ar = Dr'*Br*Dr;
	As = Ds'*Bs*Ds;
	
	Br = ABu(Rx,Rx,Br);
	Bs = ABu(Ry,Ry,Bs);
	Ar = ABu(Rx,Rx,Ar);
	As = ABu(Ry,Ry,As);
	
	[Sr,Lr] = eig(Ar,Br); Sri = inv(Sr);
	[Ss,Ls] = eig(As,Bs); Ssi = inv(Ss);
	Lfdm = diag(Lr) + diag(Ls)';
	Bm1i = 1 ./ Bm1;

elseif(slv==2) % exact solve setup

	% operators
	R  = sparse(kron(Ry,Rx));
	J  = sparse(kron(Js1d,Jr1d));
	B  = sparse(diag(reshape(Bm1,[nx1*ny1,1])));
	Bd = sparse(diag(reshape(Bmd,[nxd*nyd,1])));

	% laplace op
	G11=sparse(diag(reshape(G11,[nxd*nyd,1]))); G11 = J'*G11*J;
	G12=sparse(diag(reshape(G12,[nxd*nyd,1]))); G12 = J'*G12*J;
	G22=sparse(diag(reshape(G22,[nxd*nyd,1]))); G22 = J'*G22*J;
	G  = [G11 G12; G12 G22];
	Dr = kron(Ism1,Drm1);
	Ds = kron(Dsm1,Irm1);
	D  = [Dr;Ds];
	A  = D'*G*D;

	% advection op
	Cx = reshape(cx,[nx1*ny1,1]);
	Cy = reshape(cy,[nx1*ny1,1]);
	RRX=sparse(diag(reshape(rxm1,[nx1*ny1,1])));
	RRY=sparse(diag(reshape(rym1,[nx1*ny1,1])));
	SSX=sparse(diag(reshape(sxm1,[nx1*ny1,1])));
	SSY=sparse(diag(reshape(sym1,[nx1*ny1,1])));
	Cxd=sparse(diag(J*Cx));
	Cyd=sparse(diag(J*Cy));
	Dx = RRX*Dr + SSX*Ds;
	Dy = RRY*Dr + SSY*Ds;
	C  = J'*Bd*(Cxd*J*Dx + Cyd*J*Dy);

end

%----------------------------------------------------------------------
% LHS

function [ulhs] =  lhs_op(vlhs)
	ulhs = nu*laplace2d(vlhs,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22);
	ulhs = ulhs + vlhs.*b(1);
end
%----------------------------------------------------------------------
% RHS

function [urhs] = rhs_op(vrhs)
	urhs = -advect2d(vrhs,cx,cy,Bmd,Irm1,Ism1,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	urhs = urhs + mass2d(Bmd,Jr1d,Js1d,f); % forcing
	urhs = urhs - lhs_op(ub);              % dirichlet BC
end
%----------------------------------------------------------------------
% time advance

mesh(xm1,ym1,u);
xlabel('$$x$$');
ylabel('$$y$$');

t  = 0;
u0 = u*0;
u1 = u0;
u2 = u0;
u3 = u0;
f1 = u0;
f2 = u0;

for it=0:nt-1

	% update time
	t = t + dt;

	if(it<=3)
		[a,b] = bdfex3(i,dt);
		Lfdmi = 1 ./ (b(1) + Lfdm);
	end;
	
	% update histories
	u3=u2; u2=u1; u1 = u;
	f3=f2; f2=f1; f1 = rhs_op(u); % current u

	% form BDF rhs
	r  = a(1)*f1 +a(2)*f2 +a(3)*f3 - (b(2)*u1+b(3)*u2+b(4)*u3);
	r  = msk .* r;

	% solve
	uh = solve(r);

	u  = uh + ub;

	if(mod(it,0.1*nt)==0);
		hold off;
		mesh(xm1,ym1,u); title(['t=',num2str(t)]);pause(0.05);
	end

end
%----------------------------------------------------------------------
% solve

function [uslv] = solve(rhs)

	if(slv==0) % CG
	
		[uslv,iter,r2] = cg(rhs,0*r,1e-8,1e3);
		['cg iter:',num2str(iter),', residual:',num2str(sqrt(r2))]
	
	elseif(slv==1) % FDM
	
		uslv = (1/nu) * fdm(rhs,Bm1i,Sr,Ss,Sri,Ssi,Rx,Ry,Lfdmi);
	
	elseif(slv==2) % direct solve
	
		rhs = reshape(rhs,[nx1*ny1,1]); % rhs

		% system
		S    = R*(nu*A + b(1)*B)*R';
		rhs  = R*rhs;
		uslv = S \ rhs;
		uslv = R'*uslv;

		uslv = reshape(uslv,[nx1,ny1]);
		rhs  = reshape(rhs ,[nx1,ny1]);
	
	end

end

%======================================================================
%
%	Conjugate Gradient
%
% 	from https://en.wikipedia.org/wiki/Conjugate_gradient_method
%
%======================================================================
%
function [x,k,rsqnew] = cg(b,x0,tol,maxiter);
	x = x0;
	r = b - lhs_op(x); % r = b - Ax
	rsqold=dot2d(r,r);
	
	if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;
	
	p=r;
	for k=1:maxiter
		Ap = lhs_op(p); % Ap = A*p
		al = rsqold / dot2d(p,Ap);
		x  = x + al*p;
		r  = r - al*Ap;
		rsqnew=dot2d(r,r); if(sqrt(rsqnew) < tol); return; end;
		be = rsqnew / rsqold;
		p  = r + be*p;
		rsqold = rsqnew;
	end
end
%----------------------------------------------------------------------
end % driver
%----------------------------------------------------------------------
