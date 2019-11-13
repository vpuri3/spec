%===============================================================================
%
%	Driver function for advection diffusion equation
%
%	du/dt + \vect{c}\dot\grad{u}  = f + nu*\del^2 u
%
%   + Dirichlet/Neumann BC
%
%===============================================================================
function driver
%
%-------------------------------------------------------------------------------
%
%	/todo
%	- embed periodicity in mask using R'*R
%	- failing as nx1/ny1 > 40
%
%-------------------------------------------------------------------------------

clf; format compact; format shorte;

nx1 = 128;
ny1 = 128;
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

%-------------------------------------------------------------------------------
% geometry

[xrm1,xrp1,xsm1,xsp1,yrm1,yrp1,ysm1,ysp1] = qtrcirc(zrm1,zsm1);
[xrm2,xrp2,xsm2,xsp2,yrm2,yrp2,ysm2,ysp2] = qtrcirc(zrm2,zsm2);
[xrmd,xrpd,xsmd,xspd,yrmd,yrpd,ysmd,yspd] = qtrcirc(zrmd,zsmd);

[xm1,ym1] = gordonhall2d(xrm1,xrp1,xsm1,xsp1,yrm1,yrp1,ysm1,ysp1,zrm1,zsm1);
[xm2,ym2] = gordonhall2d(xrm2,xrp2,xsm2,xsp2,yrm2,yrp2,ysm2,ysp2,zrm2,zsm2);
[xmd,ymd] = gordonhall2d(xrmd,xrpd,xsmd,xspd,yrmd,yrpd,ysmd,yspd,zrmd,zsmd);

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xm2,ym2] = ndgrid(zrm2,zsm2);
[xmd,ymd] = ndgrid(zrmd,zsmd);

%-------------------------------------------------------------------------------
% lid driven cavity

slv=1;                           % solver --> 0: CG, 1: FDM, 2: direct

nu = 1/2e1;
vx = 0*xm1; vx(:,end) = 1;
vy = 0*xm1;
p  = 0*xm2;
fx = 0*xm1; % forcing
fy = 0*xm1;

% BC --> smooth functinos
vxb = vx;
vyb = vy;

Rxvx = Irm1(2:end-1,:);                  % dir-dir
Ryvx = Ism1(2:end-1,:);                  % dir-dir
Rxvy = Irm1(2:end-1,:);                  % dir-dir
Ryvy = Ism1(2:end-1,:);                  % dir-dir

ifconv = 1;
ifpres = 1;

T   = 10; % T=0 ==> steady
CFL = 0.5;

%------------------------------------------------------------------------------
% kovazny flow
ifkov = 1;
if(ifkov);

a = -0.5; lx = 2.5; ly=2.0;
xx = a + lx/2 * (zrm1+1) ; yy = a + ly/2 * (zsm1+1); %linear mapping
[xm1,ym1] = ndgrid(xx,yy);
xx = a + lx/2 * (zrm2+1) ; yy = a + ly/2 * (zsm2+1);
[xm2,ym2] = ndgrid(xx,yy);
xx = a + lx/2 * (zrmd+1) ; yy = a + ly/2 * (zsmd+1);
[xmd,ymd] = ndgrid(xx,yy);  

Re = 40;
Re = 10;
nu = 1/Re;

% boundary conditions
vx = 0*xm1;
vy = 0*xm1;

%left face: x=-0.5
[ue,ve] = kov_ex(xm1(1,:),ym1(1,:),Re);
vx(1,:) = ue;
vy(1,:) = ve;

%right face: x=2.0
[ue,ve] = kov_ex(xm1(end,:),ym1(end,:),Re);
vx(end,:) = ue;
vy(end,:) = ve;

%bottom face: y=-0.5
[ue,ve] = kov_ex(xm1(:,1),ym1(:,1),Re);
vx(:,1) = ue;
vy(:,1) = ve;

%top face: y=1.5
[ue,ve] = kov_ex(xm1(:,end),ym1(:,end),Re);
vx(:,end) = ue;
vy(:,end) = ve;

%p  = 0*xm2; no pressure BCs in this case
fx = 0*xm1; % forcing
fy = 0*xm1;

% BC --> smooth functinos
vxb = vx;
vyb = vy;
end
%------------------------------------------------------------------------------
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
Bm2i= 1 ./ Bm2;

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';

% laplace operator setup
g11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
g12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
g22 = Bmd .* (sxmd.*sxmd + symd.*symd);

if(slv==1) % fast diagonalization setup

	% Velocity
	Lx = max(max(xm1))-min(min(xm1));
	Ly = max(max(ym1))-min(min(ym1));
	
	Brv = (Lx/2)*diag(wrm1);
	Bsv = (Ly/2)*diag(wsm1);
	Drv = (2/Lx)*Drm1;
	Dsv = (2/Ly)*Dsm1;
	Arv = Drv'*Brv*Drv;
	Asv = Dsv'*Bsv*Dsv;
	
	Brvx = Rxvx*Brv*Rxvx';
	Bsvx = Ryvx*Bsv*Ryvx';
	Arvx = Rxvx*Arv*Rxvx';
	Asvx = Ryvx*Asv*Ryvx';
	
	Brvy = Rxvy*Brv*Rxvy';
	Bsvy = Ryvy*Bsv*Ryvy';
	Arvy = Rxvy*Arv*Rxvy';
	Asvy = Ryvy*Asv*Ryvy';
	
	[Srvx,Lrvx] = eig(Arvx,Brvx); Srivx = inv(Srvx);
	[Ssvx,Lsvx] = eig(Asvx,Bsvx); Ssivx = inv(Ssvx);
	Lvx = nu * (diag(Lrvx) + diag(Lsvx)');
	
	[Srvy,Lrvy] = eig(Arvy,Brvy); Srivy = inv(Srvy);
	[Ssvy,Lsvy] = eig(Asvy,Bsvy); Ssivy = inv(Ssvy);
	Lvy = nu * (diag(Lrvy) + diag(Lsvy)');
	
	% Pressure
	Brp  = diag(wrm2);
	Bsp  = diag(wsm2);
	Brvi = diag(1./wrm1);
	Bsvi = diag(1./wsm1);
	
	Brp = Brp*Jr12*(Drv*Brvi*Drv')*Jr12'*Brp;
	Arp = Brp*Jr12*(    Brvi     )*Jr12'*Brp;
	
	Bsp = Bsp*Js12*(Dsv*Bsvi*Dsv')*Js12'*Bsp;
	Asp = Bsp*Js12*(    Bsvi     )*Js12'*Bsp;

	[Srp,Lrp] = eig(Arp,Brp); Srip = inv(Srp);
	[Ssp,Lsp] = eig(Asp,Bsp); Ssip = inv(Ssp);
	Lp = diag(Lrp) + diag(Lsp)';
	
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

%------------------------------------------------------------------------------
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

quiver(xm1,ym1,vx,vy);grid on;
[num2str(dot(vx,Bm1.*vx)),', ',num2str(dot(vy,Bm1.*vy))]
title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);pause(0.01)

for it=1:nt

	% update histories
	t3=t2; t2=t1; t1 = t;

	p0=p;

	vx3=vx2; vx2=vx1; vx1 = vx;
	vy3=vy2; vy2=vy1; vy1 = vy;

	fvx3=fvx2; fvx2=fvx1;
	fvy3=fvy2; fvy2=fvy1;

	% pressure forcing
	[px,py] = vgradp(p0,mskvx,mskvy,Bm2,Jr12,Js12,Irm1,Ism1,Drm1...
				             ,Dsm1,rxm1,rym1,sxm1,sym1);

 	fvx1 = bdf_expl(vx,vxb,mskvx,fx,vx,vy) + px;
 	fvy1 = bdf_expl(vy,vyb,mskvy,fy,vx,vy) + py;

	t = t + dt;

	if(it<=3)
		[a,b] = bdfext3([t t1 t2 t3]);
		if(T  ==0) a=0*a; b=0*b; a(1)=1;   end; % steady
		if(slv==1) Lvxi = 1    ./ (b(1) + Lvx); % FDM
		           Lvyi = 1    ./ (b(1) + Lvy);
		           Lpi  = b(1) ./         Lp  ; end;
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
	%if(it<100) [dot(vx,Bm1.*vx),dot(vy,Bm1.*vy)],end
	if(mod(it,100)==0)
	  	quiver(xm1,ym1,vx,vy); grid on;
		%surf(xm1,ym1,vx); view(2); shading interp; title('$$v_x$$');colorbar;
	   	title(['t=',num2str(t),', Step',num2str(it),' CFL=',num2str(CFL)]);
		[dot(vx,Bm1.*vx),dot(vy,Bm1.*vy)]
		pause(0.01)
	end

end
%-------------------------------------------------------------------------------
% post process

surf(xm1,ym1,vx); view(2); shading interp; title('$$v_x$$');colorbar;
xlabel('x');
ylabel('y');
pause;
surf(xm1,ym1,vy); shading interp; view(2); title('$$v_y$$');colorbar;

%===============================================================================
%
%	Helper functions
%
%===============================================================================
% BDF - implicit OP
function [Hu] =  hmhltz(uhm,mskhm)
	Hu = nu * lapl(uhm,mskhm,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
	Hu = Hu + b(1)*mass(uhm,mskhm,Bmd,Jr1d,Js1d);
end

%-------------------------------------------------------------------------------
% BDF - explicit OP
function [Fu] = bdf_expl(uexp,ubexp,mskexp,f,cx,cy)
	if(ifconv);
		Fu = -advect(uexp,mskexp,cx,cy,Bmd,Irm1,Ism1,Jr1d... % convection
				,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	else
		Fu = 0*xm1;
	end

	Fu = Fu + mass(f,mskexp,Bmd,Jr1d,Js1d);              % forcing
	Fu = Fu - hmhltz(ubexp,1+0*mskexp);                  % dirichlet BC
end

%-------------------------------------------------------------------------------
% viscous solve
function [ux,uy] = visc_slv(rhsvx,rhsvy)

	if(slv==0) % CG
	
		ux = cg_visc(rhsvx,mskvx,0*rhsvx,1e-8,1e3);
		uy = cg_visc(rhsvy,mskvx,0*rhsvy,1e-8,1e3);
	
	elseif(slv==1) % FDM
	
		ux = fdm(rhsvx,Bm1i,Srvx,Ssvx,Srivx,Ssivx,Rxvx,Ryvx,Lvxi);
		uy = fdm(rhsvy,Bm1i,Srvy,Ssvy,Srivy,Ssivy,Rxvy,Ryvy,Lvyi);
	
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

%-------------------------------------------------------------------------------
% pressure project

function [ux,uy,p] = pres_proj(cx,cy)

	g = -qdivu(cx,cy,mskvx,mskvy,Bm2,Jr12,Js12,Irm1...
			  ,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

	if(slv==1) % FDM
		p = fdm(g,Bm2i,Srp,Ssp,Srip,Ssip,Irm2,Ism2,Lpi);
	end

	[px,py] = vgradp(p,mskvx,mskvy,Bm2,Jr12,Js12...
	                ,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

	ux = cx + 1/(b(1))*Bm1i.*px;
	uy = cy + 1/(b(1))*Bm1i.*py;

	p = p + p0;

end
%-------------------------------------------------------------------------------
% Conjugate Gradient
%
% ref https://en.wikipedia.org/wiki/Conjugate_gradient_method
%-------------------------------------------------------------------------------

function cg_pres()

end
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
%===============================================================================
end % driver
