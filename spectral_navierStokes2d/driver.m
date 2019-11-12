%======================================================================
%
%	Driver function for advection diffusion equation
%
%	du/dt + \vect{c}\dot\grad{u}  = f + nu*\del^2 u
%
%   + Dirichlet/Neumann BC
%
%======================================================================
function driver
%
%----------------------------------------------------------------------
%
%	/todo
%	- adjust mask to account for periodic BC
%	- uzawa splitting
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
[zrmd,wrmd] = zwgl (nxd-1);
[zsmd,wsmd] = zwgl (nyd-1);

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

[xrm1,xrp1,xsm1,xsp1,yrm1,yrp1,ysm1,ysp1] = qtrcirc(zrm1,zsm1);
[xrm2,xrp2,xsm2,xsp2,yrm2,yrp2,ysm2,ysp2] = qtrcirc(zrm2,zsm2);
[xrmd,xrpd,xsmd,xspd,yrmd,yrpd,ysmd,yspd] = qtrcirc(zrmd,zsmd);

[xm1,ym1] = gordonhall2d(xrm1,xrp1,xsm1,xsp1,yrm1,yrp1,ysm1,ysp1,zrm1,zsm1);
[xm2,ym2] = gordonhall2d(xrm2,xrp2,xsm2,xsp2,yrm2,yrp2,ysm2,ysp2,zrm2,zsm2);
[xmd,ymd] = gordonhall2d(xrmd,xrpd,xsmd,xspd,yrmd,yrpd,ysmd,yspd,zrmd,zsmd);

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xm2,ym2] = ndgrid(zrm2,zsm2);
[xmd,ymd] = ndgrid(zrmd,zsmd);

%----------------------------------------------------------------------
% data

slv=1;                           % solver --> 0: CG, 1: FDM, 2: direct

nu= 1e-4;
vx= (xm1-0.5).^2 + (ym1-0).^2; vx = exp(-vx/0.016);
vy= (xm1-0.5).^2 + (ym1-0).^2; vy = exp(-vy/0.016);
p = 0*xm2;
f = 0*xm1;

% BC --> smooth functinos
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
[Jm2,Jim2,rxm2,rym2,sxm2,sym2] = jac2d(xm2,ym2,Irm2,Ism2,Drm2,Dsm2);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% mass
Bm1 = Jm1.*(wrm1*wsm1');
Bm2 = Jm2.*(wrm2*wsm2');
Bmd = Jmd.*(wrmd*wsmd');
Bm1i= 1 ./ Bm1;

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';

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

Brvx = Rxvx*Br*Rxvx';
Bsvx = Ryvx*Bs*Ryvx';
Arvx = Rxvx*Ar*Rxvx';
Asvx = Ryvx*As*Ryvx';

Brvy = Rxvy*Br*Rxvy';
Bsvy = Ryvy*Bs*Ryvy';
Arvy = Rxvy*Ar*Rxvy';
Asvy = Ryvy*As*Ryvy';

[Srvx,Lrvx] = eig(Arvx,Brvx); Srivx = inv(Srvx);
[Ssvx,Lsvx] = eig(Asvx,Bsvx); Ssivx = inv(Ssvx);
Lvx = nu * (diag(Lrvx) + diag(Lsvx)');

[Srvy,Lrvy] = eig(Arvy,Brvy); Srivy = inv(Srvy);
[Ssvy,Lsvy] = eig(Asvy,Bsvy); Ssivy = inv(Ssvy);
Lvy = nu * (diag(Lrvy) + diag(Lsvy)');

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
fvx1 = vx0;
fvx2 = vx0;
% vy
vy0 = vx0;
vy1 = vx0;
vy2 = vx0;
fvy1 = vx0;
fvy2 = vx0;

quiver(xm1,ym1,vx,vy);pause(0.05);
title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);pause(0.05)

for it=1:nt

	% update histories
	t3=t2; t2=t1; t1 = t;

	p0=p;

	vx3=vx2; vx2=vx1; vx1 = vx;
	vy3=vy2; vy2=vy1; vy1 = vy;

	fvx3=fvx2; fvx2=fvx1;
	fvy3=fvy2; fvy2=fvy1;

	% pressure forcing
	[px,py] = vgradp(p0,mskvx,mskvy,Jr12,Js12,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

 	fvx1 = bdf_expl(vx,vxb,mskvx,vx,vy) + px;
 	fvy1 = bdf_expl(vy,vxb,mskvy,vx,vy) + py;

	t = t + dt;

	if(it<=3)
		[a,b] = bdfext3([t t1 t2 t3]);
		if(T  ==0) a(1)=1; b=0*b;            end; % steady
		if(slv==1) Lvxi = 1 ./ (b(1) + Lvx);      % FDM
		           Lvyi = 1 ./ (b(1) + Lvy); end;
	end;
	
	% BDF rhs
	rvx = a(1)*fvx1 +a(2)*fvx2 +a(3)*fvx3 - Bm1.*(b(2)*vx1+b(3)*vx2+b(4)*vx3);
	rvy = a(1)*fvy1 +a(2)*fvy2 +a(3)*fvy3 - Bm1.*(b(2)*vy1+b(3)*vy2+b(4)*vy3);

	% viscous solve
	[vxh,vyh] = visc_slv(rvx,rvy);

	vx  = vxh + vxb;
	vy  = vyh + vyb;

	% pressure project
	[vxh,vyh,p] = pres_proj(vx,vy);

	% vis
	if(mod(it,floor(0.1*nt))==0)
		quiver(xm1,ym1,vx,vy);
		title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);
		pause(0.05)
	end

	% keep track of norms, etc

end
%----------------------------------------------------------------------
% post process

%======================================================================
%
%	Helper functions
%
%======================================================================

% BDF - implicit OP
function [Hu] =  hmhltz(u,msk)
	Hu = nu*laplace(u,msk,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
	Hu = Hu + b(1)*mass(u,msk,Bmd,Jr1d,Js1d);
end

% BDF - explicit OP
function [Fu] = bdf_expl(u,ub,msk,cx,cy)
	Fu = -advect(u,msk,cx,cy,Bmd,Irm1,Ism1,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	Fu = Fu + mass(Bmd,Jr1d,Js1d,f); % forcing
	Fu = Fu - hmhltz(ub);            % dirichlet BC
end

% viscous solve
function [ux,uy] = visc_slv(rhsvx,rhsvy)

	if(slv==0) % CG
	
		ux = cg_visc(rhsvx,mskvx,0*rhsvx,1e-8,1e3);
		uy = cg_visc(rhsvy,mskvx,0*rhsvy,1e-8,1e3);
	
	elseif(slv==1) % FDM
	
		ux = fdm(rhsvx,Bm1i,Srvx,Ssvx,Srivx,Ssivx,Rxvx,Ryvx,Lvxi);
		uy = fdm(rhsvy,Bm1i,Srvx,Ssvx,Srivx,Ssivx,Rxvx,Ryvx,Lvyi);
	
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

% Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method

function [x,k,rsqnew] = cg_visc(b,msk,x0,tol,maxiter);
	x = x0;
	r = b - hmhltz(x,msk); % r = b - Ax
	rsqold=dot(r,r);
	
	if(sqrt(rsqold) < tol); rsqnew=rsqold; return; end;
	
	p=r;
	for k=1:maxiter
		Ap = hmhltz(p,msk); % Ap = A*p
		al = rsqold / dot(p,Ap);
		x  = x + al*p;
		r  = r - al*Ap;
		rsqnew=dot(r,r); if(sqrt(rsqnew) < tol); return; end;
		be = rsqnew / rsqold;
		p  = r + be*p;
		rsqold = rsqnew;
	end
	['cg iter:',num2str(k),', residual:',num2str(sqrt(rsqnew))]
end
%----------------------------------------------------------------------
% pressure project

function [cx,cy,p] = pres_proj(ux,uy,q0)

	g = qdivu(ux,uy,mskvx,mskvy,Jr12,Js12,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1)

	% solve for pressure
	[qx,qy] = vgradp(q,mskvx,mskvy,jr12,js12,irm1,ism1,drm1,dsm1,rxm1,rym1,sxm1,sym1)


	% adjust u 
	
end
function cg_pres()

end
%----------------------------------------------------------------------
end % driver
%----------------------------------------------------------------------
