%======================================================================
%
%	Driver function for advection diffusion equation
%
%	du/dt + \vect{c}\dot\grad{u}  = f + nu*\del^2 u
%
%   + Dirichlet/Neumann/Periodic BC
%
%======================================================================
function driver
%
%----------------------------------------------------------------------
%
%	/todo
%	- adjust mask to account for periodic BC
%
%----------------------------------------------------------------------

clf; format compact; format shorte;

nx1 = 32;
ny1 = 32;
nx2 = nx1 - 2;
ny2 = ny1 - 2;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);

[zrm1,wrm1] = zwgll(nx1-1);
[zsm1,wsm1] = zwgll(ny1-1);
[zrm2,wrm2] = zwgll(nx2-1);
[zsm2,wsm2] = zwgll(ny2-1);
[zrmd,wrmd] = zwgll(nxd-1);
[zsmd,wsmd] = zwgll(nyd-1);

Drm1 = dhat(zrm1);
Dsm1 = dhat(zsm1);
Drm2 = dhat(zrm2);
Dsm2 = dhat(zsm2);
Drmd = dhat(zrmd);
Dsmd = dhat(zsmd);

Irm1 = eye(nx1);
Ism1 = eye(ny1);
Irm2 = eye(nx2);
Ism2 = eye(ny2);
Irmd = eye(nxd);
Ismd = eye(nyd);

Jr12 = interp_mat(zrm2,zrm1); % nx1 to nx2
Js12 = interp_mat(zsm2,zsm1);
Jr1d = interp_mat(zrmd,zrm1); % nx1 to nxd
Js1d = interp_mat(zsmd,zsm1);

%----------------------------------------------------------------------
% deform geometry

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrm1,zsm1);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrm2,zsm2);
[xm2,ym2] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm2,zsm2);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrmd,zsmd);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xm2,ym2] = ndgrid(zrm2,zsm2);
[xmd,ymd] = ndgrid(zrmd,zsmd);

%----------------------------------------------------------------------
% data

slv=1;                           % solver --> 0: CG, 1: FDM, 2: direct

nu= 1e-4;
vx = (xm1-0.5).^2 + (ym1-0).^2; u = exp(-u/0.016);
vy = (xm1-0.5).^2 + (ym1-0).^2; u = exp(-u/0.016);
f = 0*xm1;

% BC
vxb = 0*xm1;
vyb = 0*xm1;

Rxvx = Irm1(2:end-1,:);                  % dir-dir
Ryvx = Ism1(2:end-1,:);                  % dir-dir
Rxvy = Irm1(2:end-1,:);                  % dir-dir
Ryvy = Ism1(2:end-1,:);                  % dir-dir

% time (steady state if T=0)
T   = 2*pi;
CFL = 0.5;

%----------------------------------------------------------------------
% setup

% time stepper
dx = min(min(diff(xm1)));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;
t  = 0;
it = 0; % time step

if(T==0); nt=1;dt=0; end; % steady

% BDF3-EXT3
a = zeros(3,1);
b = zeros(4,1);

% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% mass
Bm1 = Jm1.*(wrm1*wsm1');
Bm2 = Jm2.*(wrm2*wsm2');
Bmd = Jmd.*(wrmd*wsmd');
Bm1i= 1 ./ Bm1;

% mask off dirichlet BC
mskx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
msky = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
vxb  = (1-mskx) .* vxb;
vyb  = (1-msky) .* vyb;

% laplace operator setup
g11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
g12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
g22 = Bmd .* (sxmd.*sxmd + symd.*symd);

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
	Lfdm = nu * (diag(Lr) + diag(Ls)');

elseif(slv==2) % exact solve setup

	% operators
	R  = sparse(kron(Ry,Rx));
	J  = sparse(kron(Js1d,Jr1d));
	B  = sparse(diag(reshape(Bm1,[nx1*ny1,1])));
	Bd = sparse(diag(reshape(Bmd,[nxd*nyd,1])));

	% laplace op
	G11=sparse(diag(reshape(g11,[nxd*nyd,1]))); G11 = J'*G11*J;
	G12=sparse(diag(reshape(g12,[nxd*nyd,1]))); G12 = J'*G12*J;
	G22=sparse(diag(reshape(g22,[nxd*nyd,1]))); G22 = J'*G22*J;
	G  = [G11 G12; G12 G22];
	Dr = kron(Ism1,Drm1);
	Ds = kron(Dsm1,Irm1);
	D  = [Dr;Ds];
	A  = D'*G*D;

end

%----------------------------------------------------------------------
% time advance

t0 = 0;
t1 = 0;
t2 = 0;
% vx
vx0 = vx*0;
vx1 = vx0;
vx2 = vx0;
fx1 = vx0;
fx2 = vx0;
% vy
vy0 = vx0;
vy1 = vx0;
vy2 = vx0;
fy1 = vx0;
fy2 = vx0;

mesh(xm1,ym1,u);
title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);pause(0.05)

for it=1:nt

	vx = mask(vx,Rxvx,Ryvx);
	vy = mask(vy,Rxvy,Ryvy);

	% update histories
	t3=t2; t2=t1; t1 = t;

	vx3=vx2; vx2=vx1; vx1 = vx;
	vy3=vy2; vy2=vy1; vy1 = vy;

	fx3=fx2; fx2=fx1; fx1 = rhs_op(vx,vx,vy);
	fy3=fy2; fy2=fy1; fy1 = rhs_op(vy,vx,vy);

	t = t + dt;

	if(it<=3)
		[a,b] = bdfext3([t t1 t2 t3]);
		if(T  ==0) a(1)=1; b=0*b;              end; % steady
		if(slv==1) Lfdmi = 1 ./ (b(1) + Lfdm); end; % FDM
	end;
	
	% BDF rhs
	rx = a(1)*fx1 +a(2)*fx2 +a(3)*fx3 - Bm1.*(b(2)*vx1+b(3)*vx2+b(4)*vx3);
	ry = a(1)*fy1 +a(2)*fy2 +a(3)*fy3 - Bm1.*(b(2)*vy1+b(3)*vy2+b(4)*vy3);

	rx = mask(vx,Rxvx,Ryvx);
	ry = mask(vy,Rxvy,Ryvy);

	% viscous solve
	vxh = visc_slv(rx);
	vyh = visc_slv(ry);

	% pressure project
	[vxh,vyh] = pres_proj(vxh,vyh);

	vx  = vxh + vxb;
	vy  = vyh + vyb;

	% vis
	if(mod(it,floor(0.1*nt))==0)
		quiver(xm1,ym1,vx,vy);
		title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);
		pause(0.05)
	end

end
%----------------------------------------------------------------------
% post process

%----------------------------------------------------------------------
% viscous solve

function [uslv] = visc_slv(rslv)

	if(slv==0) % CG
	
		uslv = cg_visc(rslv,0*rslv,1e-8,1e3);
	
	elseif(slv==1) % FDM
	
		uslv = fdm(rslv,Bm1i,Sr,Ss,Sri,Ssi,Rx,Ry,Lfdmi);
	
	elseif(slv==2) % direct solve
	
		rslv = reshape(rslv,[nx1*ny1,1]);

		S    = R *(b(1)*B + nu*A)*R';
		rslv = R *rslv;
		uslv = S \rslv;
		uslv = R'*uslv;
		rslv = R'*rslv;
		
		uslv = reshape(uslv,[nx1,ny1]);
		rslv = reshape(rslv,[nx1,ny1]);
	
	end

end

%----------------------------------
% BDF - implicit OP
function [ulhs] =  lhs_op(vlhs)
	ulhs = nu*laplace2d(vlhs,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
	ulhs = ulhs + b(1)*mass2d(Bmd,Jr1d,Js1d,vlhs);

	ulhs = mask(ulhs,Rx,Ry);
end
%----------------------------------
% BDF - explicit OP
function [urhs] = rhs_op(vrhs,ux,uy)
	urhs = -advect2d(vrhs,ux,uy,Bmd,Irm1,Ism1,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	urhs = urhs + mass2d(Bmd,Jr1d,Js1d,f); % forcing
	urhs = urhs - lhs_op(ub);              % dirichlet BC

	urhs = mask(urhs,Rx,Ry);
end
%----------------------------------
% Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method

function [x,k,rsqnew] = cg_visc(b,x0,tol,maxiter);
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
	['cg iter:',num2str(k),', residual:',num2str(sqrt(rsqnew))]
end
%----------------------------------------------------------------------
% pressure project

function pres_proj(ux,uy)
end
%----------------------------------
function [v] D(ux,uy)
	[uxdx,uxdy] = grad2d(ux,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	[uydx,uydy] = grad2d(uy,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	
	v = Bm2 .* (uxdx + uydy);
end
%----------------------------------
function [ux,uy] = Dt(q)
	Bq = mass2d(Bm2,Irm2,Ism2,q);
	
end
%----------------------------------
function cg_pres()

end
%----------------------------------------------------------------------
end % driver
%----------------------------------------------------------------------
