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
%	- spectral --> spectral element
%	- Preconditioners
%		- Hlmhltz -> add viscous information
%		- Pres    -> Schwarz? Soln. to lapl eqn
%	- artificial viscosity
%
%-------------------------------------------------------------------------------
% notes
%
% field variables are defined only on the global mesh
%
%-------------------------------------------------------------------------------
clear;
clf; fig=gcf;
format compact; format shorte;

%------------
Ex  = 1;
Ey  = 1;
nx1 = 32;
ny1 = nx1;

ifvel  = 1;    % evolve velocity field per Navier-Stokes
ifpres = 0;    % project velocity field onto a div-free subspace
ifps   = 0;    % evolve passive scalar per advection diffusion eqn
%------------

nx2 = nx1 - 2;
ny2 = ny1 - 2;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);
nxp = 3*nx1;
nyp = 3*ny1;

[zrm1,wrm1] = zwgll(nx1-1); [zsm1,wsm1] = zwgll(ny1-1); % vel, scalar
[zrm2,wrm2] = zwgll(nx2-1); [zsm2,wsm2] = zwgll(ny2-1); % pres
[zrmd,wrmd] = zwgll(nxd-1); [zsmd,wsmd] = zwgll(nyd-1); % dealias
[zrmp,~   ] = zwuni(nxp-1); [zsmp,~   ] = zwgll(nyp-1); % plt

% intra-element operators
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

% global operators
Ixm1 = kron(speye(Ex),Irm1); Iym1 = kron(speye(Ey),Ism1);
Ixm2 = kron(speye(Ex),Irm2); Iym2 = kron(speye(Ey),Ism2);
Ixmd = kron(speye(Ex),Irmd); Iymd = kron(speye(Ey),Ismd);

Dxm1 = kron(speye(Ex),Drm1); Dym1 = kron(speye(Ey),Dsm1);
Dxm2 = kron(speye(Ex),Drm2); Dym2 = kron(speye(Ey),Dsm2);
Dxmd = kron(speye(Ex),Drmd); Dymd = kron(speye(Ey),Dsmd);

wxm1 = kron(ones(Ex,1),wrm1); wym1 = kron(ones(Ey,1),wsm1);
wxm2 = kron(ones(Ex,1),wrm2); wym2 = kron(ones(Ey,1),wsm2);
wxmd = kron(ones(Ex,1),wrmd); wymd = kron(ones(Ey,1),wsmd);

Jx1d = kron(speye(Ex),Jr1d); Jy1d = kron(speye(Ey),Js1d);
Jx21 = kron(speye(Ex),Jr21); Jy21 = kron(speye(Ey),Js21);
Jx1p = kron(speye(Ex),Jr1p); Jy1p = kron(speye(Ey),Js1p);
Jx2p = kron(speye(Ex),Jr2p); Jy2p = kron(speye(Ey),Js2p);

% bc
ifxprdc = 0;
ifyprdc = 0;

Rxvx = Ixm1(2:end-1,:); Ryvx = Ixm1(2:end-1,:); % restriction
Rxvy = Ixm1(2:end-1,:); Ryvy = Ixm1(2:end-1,:);
Rxps = Ixm1(2:end-1,:); Ryps = Ixm1(2:end-1,:);

Mvx   = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)'; % mask
Mvy   = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
Mps   = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';
Mbdry = 1+0*Mvx; % mask for transmiting boundary data

% global -> local operator
Qx1 = semq(Ex,nx1-1,ifxprdc); Qy1 = semq(Ey,ny1-1,ifyprdc);
Qx2 = semq(Ex,nx2-1,ifxprdc); Qy2 = semq(Ey,ny2-1,ifyprdc);
Qxd = semq(Ex,nxd-1,ifxprdc); Qyd = semq(Ey,nyd-1,ifyprdc);

%-------------------------------------------------------------------------------
% grid

% global
xm1g = semmesh(Ex,nx1,0); ym1g = semmesh(Ey,ny1,0);
xm2g = semmesh(Ex,nx2,0); ym2g = semmesh(Ey,ny2,0);
xmdg = semmesh(Ex,nxd,0); ymdg = semmesh(Ey,nyd,0);
xmpg = semmesh(Ex,nxp,0); ympg = semmesh(Ey,nyp,0);

%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xm1g,ym1g);
%[xm1g,ym1g] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xm1g,ym1g);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xm2g,ym2g);
%[xm2g,ym2g] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xm2g,ym2g);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xmdg,ymdg);
%[xmdg,ymdg] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xmdg,ymdg);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xmpg,ympg);
%[xmpg,ympg] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xmpg,ympg);

a = -0.5; lx = 2.5; ly=2.0;
xx=a+lx/2*(xm1g+1); yy=a+ly/2*(ym1g+1); [xm1g,ym1g] = ndgrid(xx,yy);
xx=a+lx/2*(xm2g+1); yy=a+ly/2*(ym2g+1); [xm2g,ym2g] = ndgrid(xx,yy);
xx=a+lx/2*(xmdg+1); yy=a+ly/2*(ymdg+1); [xmdg,ymdg] = ndgrid(xx,yy); 
xx=a+lx/2*(xmpg+1); yy=a+ly/2*(ympg+1); [xmpg,ympg] = ndgrid(xx,yy);  

% local
Qx1m = semq(Ex,nx1-1,0); Qy1m = semq(Ey,ny1-1,0);
Qx2m = semq(Ex,nx2-1,0); Qy2m = semq(Ey,ny2-1,0);
Qxdm = semq(Ex,nxd-1,0); Qydm = semq(Ey,nyd-1,0);

xm1 = ABu(Qy1m,Qx1m,xm1g); ym1 = ABu(Qy1m,Qx1m,ym1g);
xm2 = ABu(Qy2m,Qx2m,xm2g); ym2 = ABu(Qy2m,Qx2m,ym2g);
xmd = ABu(Qydm,Qxdm,xmdg); ymd = ABu(Qydm,Qxdm,ymdg);

clear xrm xrp xsm xsp yrm yrp ysm ysp Qx1m Qy1m Qx2m Qy2m Qxdm Qydm;

%-------------------------------------------------------------------------------
% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Dxm1,Dym1);
[Jm2,Jim2,rxm2,rym2,sxm2,sym2] = jac2d(xm2,ym2,Dxm2,Dym2);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Dxmd,Dymd);

% mass
Bm1  = Jm1 .* (wxm1*wym1');
Bm2  = Jm2 .* (wxm2*wym2');
Bmd  = Jmd .* (wxmd*wymd');
Bim1 = 1   ./ Bm1;

vol = dot(Bm1,1+0*Bm1);

% laplace operator setup
g11 = Bm1 .* (rxm1 .* rxm1 + rym1 .* rym1);
g12 = Bm1 .* (rxm1 .* sxm1 + rym1 .* sym1);
g22 = Bm1 .* (sxm1 .* sxm1 + sym1 .* sym1);

%-------------------------------------------------------------------------------
% kovazny

casename = 'Kovasznay Flow'; cname = 'kov';

% viscosity (velocity, passive scalar)
Re = 40;
visc0 = 1/Re;
visc1 = 0e-0;

% initial condition
vx = 0*xm1g;
vy = 0*xm1g;
ps = 0*xm1g;
pr = 0*xm2g;

% exact solution
[vxe,vye] = kov_ex(xm1g,ym1g,Re);

fvx = 0*xm1g; fvy = 0*xm1g; fps = 0*xm1g;

vxb = vxe; vyb = vye; psb = ps;

% T=0 ==> steady
T   = 20.0;
CFL = 0.5;

%------------------------------------------------------------------------------
% time stepper

dx = min(min(abs(diff(xm1))));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0) nt=1; dt=0; % steady diffusion equation
else     mov=[];     % movie
end

%------------------------------------------------------------------------------
% time advance

time = 0; 

% initialize histories
time1 = 0; time2 = 0; time3=0;
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
	%vxb = vx; vyb = vy; psb = psb;
	%fvx = fvx; fvy = fvy; fps = fps;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(it==0) ap=[1 0]; elseif(it==1) ap=[2 -1]; end
		if(T==0) a=0*a; b=0*b; a(1)=1; end; % steady heat eqn solve
	end

	% solve
	if(ifps)
		gps1 = mass(fps,Bm1,Mbdry,Qx1,Qy1) - advect(ps1,vx1,vy1,Qx1,Qy1,Bmd...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bps =       a(1)*gps1+a(2)*gps2+a(3)*gps3;
		bps = bps - mass((b(2)*ps1+b(3)*ps2+b(4)*ps3),Bm1);
		bps = bps - hlmhltz(psb,visc1,b(1),Qx1,Qy1,Bm1,Drm1,Dsm1,g11,g12,g22);
		bps = ABu(Ryps,Rxps,bps);

		psh = visc_slv(bps...
			      ,Bim1,Rxps,Ryps,visc1,b(1),Bm1,Drm1,Dsm1,g11,g12,g22);

		ps  = ABu(Ryps',Rxps',psh) + psb;
	end

	if(ifvel)
		
		gvx1 = mass(fvx,Bm1,Qx1,Qy1) - advect(vx1,vx1,vy1,Bmd...
			               		    ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
		gvy1 = mass(fvy,Bm1,Qx1,Qy1) - advect(vy1,vx1,vy1,Bmd...
						  			,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		% pressure forcing
		pr = ap(1)*pr1 + ap(2)*pr2;
		[px,py]=vgradp(pr,Bm1,Jr21,Js21,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		% viscous solve
		bvx =       a(1)*gvx1+a(2)*gvx2+a(3)*gvx3;
		bvx = bvx - mass((b(2)*vx1+b(3)*vx2+b(4)*vx3),Qx1,Qy1,Bm1);
		bvx = bvx - hlmhltz(vxb,visc0,b(1),Qx1,Qy1,Bm1,Drm1,Dsm1,g11,g12,g22);
		bvx = bvx + px;
		bvx = ABu(Ryvx,Rxvx,bvx);

		bvy =       a(1)*gvy1+a(2)*gvy2+a(3)*gvy3;
		bvy = bvy - mass((b(2)*vy1+b(3)*vy2+b(4)*vy3),Qx1,Qy1,Bm1);
		bvy = bvy - hlmhltz(vyb,visc0,b(1),Qx1,Qy1,Bm1,Drm1,Dsm1,g11,g12,g22);
		bvy = bvy + py;
		bvy = ABu(Ryvy,Rxvy,bvy);

		vxh = visc_slv(bvx...
				  ,Bim1,Rxvy,Ryvy,visc0,b(1),Bm1,Drm1,Dsm1,g11,g12,g22);
		vyh = visc_slv(bvy...
				  ,Bim1,Rxvx,Ryvx,visc0,b(1),Bm1,Drm1,Dsm1,g11,g12,g22);

		vx  = ABu(Ryvx',Rxvx',vxh) + vxb;
		vy  = ABu(Ryvy',Rxvy',vyh) + vyb;

		% pressure projection
		if(ifpres)
 			[vx,vy,pr] = pres_proj(vx,vy,pr,b(1),Bim1,Rxvx,Ryvx,Rxvy,Ryvy...
						,Bm2,Jr21,Js21,Drm1,Dsm1,rxm1,rym1,sxm1,sym1...
						,Bm1);
		end
	end

	% chk
	%---------------------------------------------------
	if(blowup(vx,vy,pr,ps,Bm1,Bm2));it, return; end;
	if(mod(it,5e1)==0 | time>=T-1e-6)

		% log
		%[it,L2(vx,Bm1),L2(vy,Bm1),L2(pr,Bm2),L2(ps,Bm1)]
		['infty kovazny normalized v-ve']
		[max(max(abs(vx-vxe))),max(max(abs(vy-vye)))] ./...
		[max(max(abs(vxe))),max(max(abs(vye)))]

		% vis
		om = vort(vx,vy,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		vxp = ABu(Js1p,Jr1p,vx);
		vyp = ABu(Js1p,Jr1p,vy);
		prp = ABu(Js2p,Jr2p,pr);
		psp = ABu(Js1p,Jr1p,ps);
		psp = ABu(Js1p,Jr1p,ps);
		omp = ABu(Js1p,Jr1p,om);

		contour(xmpg,ympg,vxp,20);
		view(2)
	   	title([casename,', t=',num2str(time,'%4.2f'),' i=',num2str(it)]);
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
