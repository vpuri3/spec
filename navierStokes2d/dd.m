%===============================================================================
%
%	Driver function for Navier Stokes equation
%
%	\partial_t u + (u\dot\grad) u = -\grad p + nu*\del^2 u + f
%   				  \grad\dot u = 0
%
%   + Dirichlet/Neumann/Periodic BC
%
%===============================================================================
%function driver
%
%-------------------------------------------------------------------------------
%
%	/todo
%	- add references (Fischer JCP 97)
%
%-------------------------------------------------------------------------------
clear;
clf; fig=gcf;
format compact; format shorte;

%------------
ifannulus=1;
ifrb     =1;

nx1 = 50;
ny1 = 50;

slv=0; % 0: CG, 1: FDM

ifvel  = 1;    % evolve velocity field per Navier-Stokes
ifpres = 1;    % project velocity field onto a div-free subspace
ifps   = 1;    % evolve passive scalar per advection diffusion
%------------

nx2 = nx1 - 2;
ny2 = ny1 - 2;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);
nxp = 10*nx1;
nyp = 10*ny1;

[zrm1,wrm1] = zwgll(nx1-1); [zsm1,wsm1] = zwgll(ny1-1); % vel, scalar
[zrm2,wrm2] = zwgll(nx2-1); [zsm2,wsm2] = zwgll(ny2-1); % pres
[zrmd,wrmd] = zwgll(nxd-1); [zsmd,wsmd] = zwgll(nyd-1); % dealias
[zrmp,~   ] = zwuni(nxd-1); [zsmp,~   ] = zwgll(nyp-1); % plt

Drm1 = dhat(zrm1); Dsm1 = dhat(zsm1);
Drm2 = dhat(zrm2); Dsm2 = dhat(zsm2);
Drmd = dhat(zrmd); Dsmd = dhat(zsmd);

Irm1 = speye(nx1); Ism1 = speye(ny1);
Irm2 = speye(nx2); Ism2 = speye(ny2);
Irmd = speye(nxd); Ismd = speye(nyd);

Jr1d = interp_mat(zrmd,zrm1); Js1d = interp_mat(zsmd,zsm1); % vel  -> dealias
Jr21 = interp_mat(zrm1,zrm2); Js21 = interp_mat(zsm1,zsm2); % pres -> vel
Jr1p = interp_mat(zrmp,zrm1); Js1p = interp_mat(zsmp,zsm1); % vel  -> plt
Jr2p = interp_mat(zrmp,zrm2); Js2p = interp_mat(zsmp,zsm2); % pres -> plt

%-------------------------------------------------------------------------------
% geometry

%[xm1,ym1] = ndgrid(zrm1,zsm1);
%[xm2,ym2] = ndgrid(zrm2,zsm2);
%[xmd,ymd] = ndgrid(zrmd,zsmd);
%[xmp,ymp] = ndgrid(zrmp,zsmp);
%
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(zrm1,zsm1);
%[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(zrm2,zsm2);
%[xm2,ym2] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm2,zsm2);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(zrmd,zsmd);
%[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(zrmp,zsmp);
%[xmp,ymp] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmp,zsmp);

kc = 3.117; a = 0; Lx = 2*pi/kc; Ly=1;
xx = a + Lx/2 * (zrm1+1); yy = a + Ly/2 * (zsm1+1); [xm1,ym1]=ndgrid(xx,yy);
xx = a + Lx/2 * (zrm2+1); yy = a + Ly/2 * (zsm2+1); [xm2,ym2]=ndgrid(xx,yy);
xx = a + Lx/2 * (zrmd+1); yy = a + Ly/2 * (zsmd+1); [xmd,ymd]=ndgrid(xx,yy);
xx = a + Lx/2 * (zrmp+1); yy = a + Ly/2 * (zsmp+1); [xmp,ymp]=ndgrid(xx,yy);

if(ifannulus)
r0=1; r1=2; span=2*pi;
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrm1,zsm1);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrm2,zsm2);
[xm2,ym2] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm2,zsm2);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrmd,zsmd);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrmp,zsmp);
[xmp,ymp] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmp,zsmp);
end

if(ifannulus)
	Lx = 0.5*span*(r1*r1-r0*r0);
	Ly = r1-r0;
else
	Lx = abs(sum(xm1(end,:)-xm1(1,:))) / nx1;
	Ly = abs(sum(ym1(:,end)-ym1(:,1))) / ny1;
end

if(0)
%------------------------------------------------------------------------------
casename = 'Annulus'; cname = 'ans';

%visc1 = 1e-1; vx=0*ym1; vy=0*xm1; fps=0*xm1; % pure diffusion
%ps=1+0*xm1; ps(:,2:end)=0;
%vxb = vx; vyb = vy; psb = ps;
%r = sqrt(xm1.*xm1 + ym1.*ym1); a = atan2(ym1,xm1);

visc0 = 1/1e2;
visc1 = 0e-0;

vx = -ym1;%vx(:,2:end) = 0;
vy =  xm1;%vy(:,2:end) = 0;
pr =0*xm2;
ps =0*xm1;

fvx = 0*xm1; fvy = 0*xm1; fps = 0*xm1;
vxb = vx; vyb = vy; psb = ps;

Rxvx = Irm1(2:end-1,:); Ryvx = Ism1(2:end-1,:);
Rxvy = Irm1(2:end-1,:); Ryvy = Ism1(2:end-1,:);
Rxps = Irm1(2:end-1,:); Ryps = Ism1(2:end-1,:);

ifxperiodic = 1;
ifyperiodic = 0;

T   = 10.0;
CFL = 0.1;
end
%------------------------------------------------------------------------------
casename = 'RB'; cname = 'rb';

Ra = 1707.76;
Pr = 0.71; %air

visc0 = Pr;
visc1 = 1;

% visc0 = 0.023833; %Ra=1e4;Pr=0.71
% visc1 = 0.033567;

% visc0 = 2.2833e-3; %Ra=1e6;Pr=0.71
% visc1 = 3.3567e-3;

% initial condition
vx  = -xm1*0;
vy  =  ym1*0;
pr  = 0*xm2;
ps  = 0*xm1; ps(:,1)=1;

fvx = 0*xm1; fvy = 0*xm1; fps = 0*xm1; 
vxb = vx; vyb = vy; psb = ps;

% Restrictions
Rxvx = Irm1(2:end-1,:); Rxvy = Irm1(2:end-1,:);
Ryvx = Ism1(2:end-1,:); Ryvy = Ism1(2:end-1,:);
Rxps = Irm1(1:end-0,:); Ryps = Ism1(2:end-1,:);

ifxperiodic = 1;
ifyperiodic = 0;

T   = 4.0;
CFL = 0.5;

%------------------------------------------------------------------------------
% setup

% periodic BC through restriction matrices
if(ifxperiodic)
	Rxvx = Irm1(1:end-1,:); Rxvx(1,end)=1; Rxvy = Rxvx; Rxps = Rxvx;
	Rxpr = Irm2(1:end-1,:); Rxpr(1,end)=1;
end;

if(ifyperiodic)
	Ryvx = Ism1(1:end-1,:); Ryvx(1,end)=1; Ryvy = Ryvx; Ryps = Ryvx;
	Rypr = Ism2(1:end-1,:); Rypr(1,end)=1;
end;

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
mskps = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';

% time stepper
dx = min(min(abs(diff(xm1))));
if(ifannulus) dx = span*r0/(nx1*nx1);end
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0) nt=1; dt=0; % steady
else     mov=[];     % movie
end

% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
[Jm2,Jim2,rxm2,rym2,sxm2,sym2] = jac2d(xm2,ym2,Irm2,Ism2,Drm2,Dsm2);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% mass
Bm1  = Jm1 .* (wrm1*wsm1');
Bm2  = Jm2 .* (wrm2*wsm2');
Bmd  = Jmd .* (wrmd*wsmd');
Bim1 = 1   ./ Bm1;

vol = dot(Bm1,1+0*Bm1);

% laplace operator setup
g11 = Bm1 .* (rxm1 .* rxm1 + rym1 .* rym1);
g12 = Bm1 .* (rxm1 .* sxm1 + rym1 .* sym1);
g22 = Bm1 .* (sxm1 .* sxm1 + sym1 .* sym1);

%------------------------------------------------------------------------------
% fast diagonalization setup

% Velocity
Bxv = (Lx/2)*diag(wrm1); Byv = (Ly/2)*diag(wsm1);
Dxv = (2/Lx)*Drm1;       Dyv = (2/Ly)*Dsm1;
Axv = Dxv'*Bxv*Dxv;      Ayv = Dyv'*Byv*Dyv;

Bxvx = Rxvx*Bxv*Rxvx'; Byvx = Ryvx*Byv*Ryvx';
Axvx = Rxvx*Axv*Rxvx'; Ayvx = Ryvx*Ayv*Ryvx';

Bxvy = Rxvy*Bxv*Rxvy'; Byvy = Ryvy*Byv*Ryvy';
Axvy = Rxvy*Axv*Rxvy'; Ayvy = Ryvy*Ayv*Ryvy';

[Sxvx,Lxvx] = eig(Axvx,Bxvx);
[Syvx,Lyvx] = eig(Ayvx,Byvx);
Sxvx=Sxvx*diag(1./sqrt(diag(Sxvx'*Bxvx*Sxvx)));
Syvx=Syvx*diag(1./sqrt(diag(Syvx'*Byvx*Syvx)));
Lvx = visc0 * (diag(Lxvx) + diag(Lyvx)');

[Sxvy,Lxvy] = eig(Axvy,Bxvy);
[Syvy,Lyvy] = eig(Ayvy,Byvy);
Sxvy=Sxvy*diag(1./sqrt(diag(Sxvy'*Bxvy*Sxvy)));
Syvy=Syvy*diag(1./sqrt(diag(Syvy'*Byvy*Syvy)));
Lvy = visc0 * (diag(Lxvy) + diag(Lyvy)');

% Passive Scalar
Bxps = Rxps*Bxv*Rxps'; Byps = Ryps*Byv*Ryps';
Axps = Rxps*Axv*Rxps'; Ayps = Ryps*Ayv*Ryps';

[Sxps,Lxps] = eig(Axps,Bxps);
[Syps,Lyps] = eig(Ayps,Byps);
Sxps=Sxps*diag(1./sqrt(diag(Sxps'*Bxps*Sxps)));
Syps=Syps*diag(1./sqrt(diag(Syps'*Byps*Syps)));
Lps = visc1 * (diag(Lxps) + diag(Lyps)');

% Pressure
Myvx = Ryvx'*Ryvx; Mxvx = Rxvx'*Rxvx;
Myvy = Ryvy'*Ryvy; Mxvy = Rxvy'*Rxvy;
Bxiv = diag(1./diag(Bxv));
Byiv = diag(1./diag(Byv));

Byp = Js21'*Byv*(    (Myvx*Byiv*Myvx)     )*Byv*Js21; % attack vx
Axp = Jr21'*Bxv*(Dxv*(Mxvx*Bxiv*Mxvx)*Dxv')*Bxv*Jr21;
Ayp = Js21'*Byv*(Dyv*(Myvy*Byiv*Myvy)*Dyv')*Byv*Js21; % attack vy
Bxp = Jr21'*Bxv*(    (Mxvy*Bxiv*Mxvy)     )*Bxv*Jr21;

[Sxpr,Lxpr] = eig(Axp,Bxp);
[Sypr,Lypr] = eig(Ayp,Byp);
Sxpr=Sxpr*diag(1./sqrt(diag(Sxpr'*Bxp*Sxpr)));
Sypr=Sypr*diag(1./sqrt(diag(Sypr'*Byp*Sypr)));
Lpr  = diag(Lxpr) + diag(Lypr)';
Lipr = 1 ./ Lpr;
Lipr(find(abs(Lipr)>1e10)) = 0;

%------------------------------------------------------------------------------
ifdebug=0;
if(ifdebug)
%J21  = kron(Js21,Jr21);
%Bv   = kron(Byv ,Bxv );
%Mvx  = kron(Myvx,Mxvx);
%Mvy  = kron(Myvy,Mxvy);
%MM   = [Mvx,zeros(nx1*ny1);zeros(nx1*ny1),Mvy];
%DDb  = J21'*Bv*[kron(Ism1,Dxv),kron(Dyv,Irm1)]; % DD_bar
%Biv  = kron(Byiv,Bxiv);
%BBiv = kron(eye(2),Biv);

%E = DDb*MM*BBiv*MM*DDb';
%F = kron(Byp,Axp)+kron(Ayp,Bxp);  ['err in forming E'],max(max(abs(F-E)))
%e=sort(eig(E)); ['e.vals of exp  mat'],e(1:6)'
%e=sort(eig(F)); ['e.vals of kron mat'],e(1:6)'

  % operators
R   = sparse(kron(Ryps,Rxps));
J1d = sparse(kron(Js1d,Jr1d));
B1  = sparse(diag(reshape(Bm1,[nx1*ny1,1])));
Bd  = sparse(diag(reshape(Bmd,[nxd*nyd,1])));
B1i = sparse(diag(reshape(1./Bm1,[nx1*ny1,1])));

 % lapl op
G11=sparse(diag(reshape(g11,[nx1*ny1,1])));
G12=sparse(diag(reshape(g12,[nx1*ny1,1])));
G22=sparse(diag(reshape(g22,[nx1*ny1,1])));
G  = [G11 G12; G12 G22];
Dr = kron(Ism1,Drm1);
Ds = kron(Dsm1,Irm1);
D  = [Dr;Ds];
A  = D'*G*D;

 % advec op
Cx = reshape(vx,[nx1*ny1,1]);
Cy = reshape(vy,[nx1*ny1,1]);
RRX=sparse(diag(reshape(rxm1,[nx1*ny1,1])));
RRY=sparse(diag(reshape(rym1,[nx1*ny1,1])));
SSX=sparse(diag(reshape(sxm1,[nx1*ny1,1])));
SSY=sparse(diag(reshape(sym1,[nx1*ny1,1])));
Cxd=sparse(diag(J1d*Cx));
Cyd=sparse(diag(J1d*Cy));
Dx = RRX*Dr + SSX*Ds;
Dy = RRY*Dr + SSY*Ds;
C  = J1d'*Bd*(Cxd*J1d*Dx+Cyd*J1d*Dy);
plot(eig(full(C)),'o'),pause
end
%------------------------------------------------------------------------------
% time advance

time = 0; 

% initialize histories
time1 = time*0; time2 = 0; time3=0;
vx1 = vx*0; vx2 = vx1; vx3 = vx2; gvx1 = vx1; gvx2 = vx1; gvx3 = vx1;
vy1 = vy*0; vy2 = vy1; vy3 = vy2; gvy1 = vy1; gvy2 = vy1; gvy3 = vy1;
ps1 = ps*0; ps2 = ps1; ps3 = ps2; gps1 = ps1; gps2 = ps1; gps3 = ps1;
pr1 = pr*0; pr2 = pr*0;

for it=1:nt

	% update histories
	time3=time2; time2=time1; time1=time; time = time + dt;

	vx3=vx2; vx2=vx1; vx1=vx; gvx3=gvx2; gvx2=gvx1;
	vy3=vy2; vy2=vy1; vy1=vy; gvy3=gvy2; gvy2=gvy1;
	ps3=ps2; ps2=ps1; ps1=ps; gps3=gps2; gps2=gps1;
			 pr2=pr1; pr1=pr;

	% update BC, forcing
	if(ifrb)
    	fvy = Ra*Pr*ps;
	end
	%vxb = vxb; vyb = vyb; psb = psb;
	%fvx = fvx; fvy = fvy; fps = fps;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(it==0) ap=[1 0]; elseif(it==1) ap=[2 -1]; end
		if(T==0) a=0*a; b=0*b; a(1)=1; end; % steady solve
		Livx = 1 ./ (b(1) + Lvx);	  	    % FDM
		Livy = 1 ./ (b(1) + Lvy);
		Lips = 1 ./ (b(1) + Lps);
	end

	% solve
	if(ifps)
		gps1 = mass(fps,Bm1,Irm1,Ism1) - advect(ps1,vx1,vy1,Bmd,Irm1,Ism1...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bps =       a(1)*gps1+a(2)*gps2+a(3)*gps3;
		bps = bps - mass((b(2)*ps1+b(3)*ps2+b(4)*ps3),Bm1,Irm1,Ism1);
		bps = bps - hlmhltz(psb,visc1,b(1),Bm1,Irm1,Ism1,Drm1,Dsm1,g11,g12,g22);
		bps = ABu(Ryps,Rxps,bps);

		psh = visc_slv(bps,Sxps,Syps,Lips,slv...
			      ,Bim1,Rxps,Ryps,visc1,b(1),Bm1,Irm1,Ism1,Drm1,Dsm1,g11,g12,g22);

		ps  = ABu(Ryps',Rxps',psh) + psb;
	end

	if(ifvel)
		
		gvx1 = mass(fvx,Bm1,Irm1,Ism1) - advect(vx1,vx1,vy1,Bmd,Irm1,Ism1...
			               		    ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
		gvy1 = mass(fvy,Bm1,Irm1,Ism1) - advect(vy1,vx1,vy1,Bmd,Irm1,Ism1...
						  			,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		% pressure forcing
		pr = ap(1)*pr1 + ap(2)*pr2;
		[px,py]=vgradp(pr,Bm1,Jr21,Js21,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		% viscous solve
		bvx =       a(1)*gvx1+a(2)*gvx2+a(3)*gvx3;
		bvx = bvx - mass((b(2)*vx1+b(3)*vx2+b(4)*vx3),Bm1,Irm1,Ism1);
		bvx = bvx - hlmhltz(vxb,visc0,b(1),Bm1,Irm1,Ism1,Drm1,Dsm1,g11,g12,g22);
		bvx = bvx + px;
		bvx = ABu(Ryvx,Rxvx,bvx);

		bvy =       a(1)*gvy1+a(2)*gvy2+a(3)*gvy3;
		bvy = bvy - mass((b(2)*vy1+b(3)*vy2+b(4)*vy3),Bm1,Irm1,Ism1);
		bvy = bvy - hlmhltz(vyb,visc0,b(1),Bm1,Irm1,Ism1,Drm1,Dsm1,g11,g12,g22);
		bvy = bvy + py;
		bvy = ABu(Ryvy,Rxvy,bvy);

		vxh = visc_slv(bvx,Sxvx,Syvx,Livx,slv...
				  ,Bim1,Rxvy,Ryvy,visc0,b(1),Bm1,Irm1,Ism1,Drm1,Dsm1,g11,g12,g22);
		vyh = visc_slv(bvy,Sxvy,Syvy,Livy,slv...
				  ,Bim1,Rxvx,Ryvx,visc0,b(1),Bm1,Irm1,Ism1,Drm1,Dsm1,g11,g12,g22);

		vx  = ABu(Ryvx',Rxvx',vxh) + vxb;
		vy  = ABu(Ryvy',Rxvy',vyh) + vyb;

		% pressure projection
		if(ifpres)
 			[vx,vy,pr] = pres_proj(vx,vy,pr,b(1),Bim1,Rxvx,Ryvx,Rxvy,Ryvy,slv...
						,Bm2,Jr21,Js21,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1...
						,Sxpr,Sypr,Lipr,Bm1);
		end
	end

	% chk
	%---------------------------------------------------

	%mesh(xm1,ym1,ps),drawnow
	if(blowup(vx,vy,pr,ps,Bm1,Bm2));it, return; end;
	if(mod(it,5e1)==0 | time>=T-1e-6)

		% log
		[it,L2(vx,Bm1),L2(vy,Bm1),L2(pr,Bm2),L2(ps,Bm1)]

		% vis
		om = vort(vx,vy,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		vxp = ABu(Js1p,Jr1p,vx);
		vyp = ABu(Js1p,Jr1p,vy);
		prp = ABu(Js2p,Jr2p,pr);
		psp = ABu(Js1p,Jr1p,ps);
		psp = ABu(Js1p,Jr1p,ps);
		omp = ABu(Js1p,Jr1p,om);

		%mesh(xmp,ymp,psp); view(3);%max(max(abs(psp)))
		contour(xmp,ymp,vyp,20);

	   	title([casename,', t=',num2str(time,'%4.2f'),' i=',num2str(it)]);
		colorbar
		drawnow
		if(T ~=0) mov = [mov,getframe(fig)]; end
	end
	%---------------------------------------------------
end
%-------------------------------------------------------------------------------
% post process

['Finished Timestepping']
%['Energy in vx,vy,pr,ps'],[L2(vx,Bm1),L2(vy,Bm1),L2(pr,Bm2),L2(ps,Bm1)]

% play movie
%movie(fig,mov,-2,40);

% save gif

%gname = [cname,'.gif'];
%fps   = 40;
%mov   = [mov,flip(mov)];
%
%for i=1:length(mov)
%	f = mov(i);
%	[img,cmap] = rgb2ind(f.cdata,256);
%	if i==1 imwrite(img,cmap,gname,'gif','DelayTime',1/fps,'LoopCount',Inf)
%	else imwrite(img,cmap,gname,'gif','WriteMode','append','DelayTime',1/fps)
%	end
%end

%===============================================================================
%end % driver
%===============================================================================
