%===============================================================================
%
%	Driver function for advection diffusion equation
%
%	du/dt + \vect{c}\dot\grad{u}  = f + nu*\del^2 u
%
%   + Dirichlet/Neumann/Periodic BC
%
%===============================================================================
function driver
%
%-------------------------------------------------------------------------------
%
%	/todo
%	- problem w pressure FDM
%	- walsch setup
%	- add references
%
%-------------------------------------------------------------------------------

clf; format compact; format shorte;

%------------
ifkov = 1; % exponential convergence achieved
ifLDC = 0; % blowing up
ifwls = 0; % /todo

nx1 = 32;
ny1 = 32;

slv=1; % 0: CG /todo, 1: FDM
%------------

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

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xm2,ym2] = ndgrid(zrm2,zsm2);
[xmd,ymd] = ndgrid(zrmd,zsmd);

%-------------------------------------------------------------------------------
% lid driven cavity
if(ifLDC)

% viscosity (velocity, passive scalar)
visc0 = 1e-2;
visc1 = 1e-2;

% initial condition
vx  = 0*xm1; vx(:,end)=1;
vy  = 0*xm1;
ps  = 0*xm1; ps(:,end)=1;
pr  = 0*xm2;

% forcing
fvx = 0*xm1;
fvy = 0*xm1;
fps = 0*xm1;%fps = sin(pi*xm1).*sin(pi*ym1); pse = fps/2/pi/pi/visc1;

% BC
vxb = vx;
vyb = vy;
psb = ps;

% Restrictions
Rxvx = Irm1(2:end-1,:); % vx             % dir-dir
Ryvx = Ism1(2:end-1,:);                  % dir-dir
Rxvy = Irm1(2:end-1,:); % vy             % dir-dir
Ryvy = Ism1(2:end-1,:);                  % dir-dir
Rxps = Irm1(2:end-1,:); % ps             % dir-dir
Ryps = Ism1(2:end-1,:);                  % dir-dir

ifxperiodic = 0;
ifyperiodic = 0;

ifvel  = 1;    % evolve velocity field
ifpres = 1;    % project velocity field onto a div-free subspace
ifps   = 0;    % evolve passive scalar per advection diffusion

% T=0 ==> steady
T   = 10;
CFL = 0.5;

end
%-------------------------------------------------------------------------------
% kovazny

if(ifkov)

a = -0.5; lx = 2.5; ly=2.0;
xx = a + lx/2 * (zrm1+1) ; yy = a + ly/2 * (zsm1+1); [xm1,ym1] = ndgrid(xx,yy);
xx = a + lx/2 * (zrm2+1) ; yy = a + ly/2 * (zsm2+1); [xm2,ym2] = ndgrid(xx,yy);
xx = a + lx/2 * (zrmd+1) ; yy = a + ly/2 * (zsmd+1); [xmd,ymd] = ndgrid(xx,yy);  

% viscosity (velocity, passive scalar)
Re = 40;
visc0 = 1/Re;
visc1 = 0e-0;

% initial condition
vx = 0*xm1;
vy = 0*xm1;
ps = 0*xm1;
pr = 0*xm2;

% exact solution
[vxe,vye] = kov_ex(xm1,ym1,Re);
%left face: x=-0.5
vx(1,:) = vxe(1,:);
vy(1,:) = vye(1,:);
%right face: x=2.0
vx(end,:) = vxe(end,:);
vy(end,:) = vye(end,:);
%bottom face: y=-0.5
vx(:,1) = vxe(:,1);
vy(:,1) = vye(:,1);
%top face: y=1.5
vx(:,end) = vxe(:,end);
vy(:,end) = vye(:,end);

% forcing
fvx = 0*xm1;
fvy = 0*xm1;
fps = 0*xm1;

% BC
vxb = vx;
vyb = vy;
psb = ps;

% Restrictions
Rxvx = Irm1(2:end-1,:); % vx             % dir-dir
Ryvx = Ism1(2:end-1,:);                  % dir-dir
Rxvy = Irm1(2:end-1,:); % vy             % dir-dir
Ryvy = Ism1(2:end-1,:);                  % dir-dir
Rxps = Irm1(2:end-1,:); % ps             % dir-dir
Ryps = Ism1(2:end-1,:);                  % dir-dir

ifxperiodic = 0;
ifyperiodic = 0;

ifvel  = 1;    % evolve velocity field
ifps   = 0;    % evolve passive scalar per advection diffusion
ifpres = 1;    % project velocity field onto a div-free subspace

% T=0 ==> steady
T   = 5.0;
CFL = 0.5;

end
%------------------------------------------------------------------------------
% Walsch
if(ifwls)

% viscosity (velocity, passive scalar)
visc0 = 1e-2;
visc1 = 1e-2;

% initial condition
vx  = 0*xm1; vx(:,end)=1;
vy  = 0*xm1;
ps  = 0*xm1; ps(:,end)=1;
pr  = 0*xm2;

% forcing
fvx = 0*xm1;
fvy = 0*xm1;
fps = 0*xm1;%fps = sin(pi*xm1).*sin(pi*ym1); pse = fps/2/pi/pi/visc1;

% BC
vxb = vx;
vyb = vy;
psb = ps;

% Restrictions
Rxvx = Irm1(2:end-1,:); % vx             % dir-dir
Ryvx = Ism1(2:end-1,:);                  % dir-dir
Rxvy = Irm1(2:end-1,:); % vy             % dir-dir
Ryvy = Ism1(2:end-1,:);                  % dir-dir
Rxps = Irm1(2:end-1,:); % ps             % dir-dir
Ryps = Ism1(2:end-1,:);                  % dir-dir

ifxperiodic = 0;
ifyperiodic = 0;

ifvel  = 1;    % evolve velocity field
ifpres = 1;    % project velocity field onto a div-free subspace
ifps   = 0;    % evolve passive scalar per advection diffusion

% T=0 ==> steady
T   = 10;
CFL = 0.5;

end
%------------------------------------------------------------------------------
% setup

% periodic BC through restriction matrices
if(ifxperiodic)
	Rxvx = [eye(nx1-1),[1;zeros(nx1-2,1)]];
	Rxvy = Rxvx;
	Rxps = Rxvx;
	Rxpr = [eye(nx2-1),[1;zeros(nx2-2,1)]];
end;

if(ifyperiodic)
	Ryvx = [eye(ny1-1),[1;zeros(ny1-2,1)]];
	Ryvy = Ryvx;
	Ryps = Ryvx;
	Rypr = [eye(ny2-1),[1;zeros(ny2-2,1)]];
end;

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
mskps = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';

% time stepper
dx = min(min(diff(xm1)));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0); nt=1;dt=0; end; % steady

% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
[Jm2,Jim2,rxm2,rym2,sxm2,sym2] = jac2d(xm2,ym2,Irm2,Ism2,Drm2,Dsm2);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% mass
Bm1  = Jm1 .* (wrm1*wsm1');
Bm2  = Jm2 .* (wrm2*wsm2');
Bmd  = Jmd .* (wrmd*wsmd');
Bim1 = 1   ./ Bm1;
Bim2 = 1   ./ Bm2;

% laplace operator setup
g11 = Bmd .* (rxmd .* rxmd + rymd .* rymd);
g12 = Bmd .* (rxmd .* sxmd + rymd .* symd);
g22 = Bmd .* (sxmd .* sxmd + symd .* symd);

%------------------------------------------------------------------------------
% fast diagonalization setup
Lx = max(max(xm1))-min(min(xm1));
Ly = max(max(ym1))-min(min(ym1));

% Velocity
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

[Srvx,Lrvx] = eig(Arvx,Brvx);
[Ssvx,Lsvx] = eig(Asvx,Bsvx);
Srvx=Srvx*diag(1./sqrt(diag(Srvx'*Brvx*Srvx)));
Ssvx=Ssvx*diag(1./sqrt(diag(Ssvx'*Bsvx*Ssvx)));
Lvx = visc0 * (diag(Lrvx) + diag(Lsvx)');

[Srvy,Lrvy] = eig(Arvy,Brvy);
[Ssvy,Lsvy] = eig(Asvy,Bsvy);
Srvy=Srvy*diag(1./sqrt(diag(Srvy'*Brvy*Srvy)));
Ssvy=Ssvy*diag(1./sqrt(diag(Ssvy'*Bsvy*Ssvy)));
Lvy = visc0 * (diag(Lrvy) + diag(Lsvy)');

% Passive Scalar

Brps = Rxps*Brv*Rxps';
Bsps = Ryps*Bsv*Ryps';
Arps = Rxps*Arv*Rxps';
Asps = Ryps*Asv*Ryps';

[Srps,Lrps] = eig(Arps,Brps);
[Ssps,Lsps] = eig(Asps,Bsps);
Srps=Srps*diag(1./sqrt(diag(Srps'*Brps*Srps)));
Ssps=Ssps*diag(1./sqrt(diag(Ssps'*Bsps*Ssps)));
Lps = visc1 * (diag(Lrps) + diag(Lsps)');

% Pressure
Jr21 = interp_mat(zrm1,zrm2); % nx2 to nx1
Js21 = interp_mat(zsm1,zsm2);

Brpr = (Lx/2)*diag(wrm2);
Bspr = (Ly/2)*diag(wsm2);
Briv = diag(1./diag(Brv));
Bsiv = diag(1./diag(Bsv));

Msvx = Ryvx'*Ryvx; % mask matrices
Mrvx = Rxvx'*Rxvx;
Msvy = Ryvy'*Ryvy;
Mrvy = Rxvy'*Rxvy;

Bsp = Bspr*Js21'*(    (Msvx*Bsiv)     )*Js21*Bspr; % attack vx
Arp = Brpr*Jr21'*(Drv*(Mrvx*Briv)*Drv')*Jr21*Brpr;
Asp = Bspr*Js21'*(Dsv*(Msvy*Bsiv)*Dsv')*Js21*Bspr; % attack vy
Brp = Brpr*Jr21'*(    (Mrvy*Briv)     )*Jr21*Brpr;
[Srpr,Lrpr] = eig(Arp,Brp);
[Sspr,Lspr] = eig(Asp,Bsp);

Bspr = Bspr*Js12*(    Bsiv     )*Js12'*Bspr; % attack vx
Arpr = Brpr*Jr12*(Drv*Briv*Drv')*Jr12'*Brpr;
Aspr = Bspr*Js12*(Dsv*Bsiv*Dsv')*Js12'*Bspr; % attack vy
Brpr = Brpr*Jr12*(    Briv     )*Jr12'*Brpr;
[Srpr,Lrpr] = eig(Arpr,Brpr);
[Sspr,Lspr] = eig(Aspr,Bspr);

Srpr=Srpr*diag(1./sqrt(diag(Srpr'*Brpr*Srpr)));
Sspr=Sspr*diag(1./sqrt(diag(Sspr'*Bspr*Sspr)));
Lpr = diag(Lrpr) + diag(Lspr)';
%diag(Lrpr)',diag(Lspr)',pause

%------------------------------------------------------------------------------
% time advance

time = 0;

% initialize histories
time0 = 0;
time1 = 0;
time2 = 0;
% vx
vx0  = vx*0;
vx1  = vx0;
vx2  = vx0;
gvx1 = vx0;
gvx2 = vx0;
% vy
vy0  = vy*0;
vy1  = vy0;
vy2  = vy0;
gvy1 = vy0;
gvy2 = vy0;
% pr
pr1  = pr*0;
% ps
ps0  = ps*0;
ps1  = ps0;
ps2  = ps0;
gps1 = ps0;
gps2 = ps0;

for it=1:nt

	time3=time2; time2=time1; time1=time;
	time = time + dt;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(T  ==0) a=0*a; b=0*b; a(1)=1;   end; % steady
		if(slv==1) Livx = 1    ./ (b(1) + Lvx); % FDM
		           Livy = 1    ./ (b(1) + Lvy);
		           Lips = 1    ./ (b(1) + Lps);
		           Lipr = b(1) ./ (       Lpr); end;
	end

	if(ifvel)
		% pressure forcing
		if(ifpres)
			pr1 = pr;
			[px,py] = vgradp(pr1,Bm2,Jr12,Js12,Irm1,Ism1...
	 						,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
		else
			px = 0*vx;
			py = 0*vy;
		end

		% vx solve
		 vx3= vx2;  vx2= vx1;  vx1 = vx;
		gvx3=gvx2; gvx2=gvx1;
		
		gvx1 = mass(fvx,Bmd,Jr1d,Js1d) - advect(vx,vx,vy,Bmd,Irm1,Ism1...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bvx =       a(1)*gvx1+a(2)*gvx2+a(3)*gvx3;
		bvx = bvx - mass((b(2)*vx1+b(3)*vx2+b(4)*vx3),Bmd,Jr1d,Js1d);
		bvx = bvx - hmhltz(vxb,visc0,b(1),Bmd,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
		bvx = bvx + px;
		bvx = ABu(Ryvx,Rxvx,bvx);

		vxh = visc_slv(bvx,Srvx,Ssvx,Livx,slv);
		vx  = ABu(Ryvx',Rxvx',vxh) + vxb;

		% vy solve
		 vy3= vy2;  vy2= vy1;  vy1 = vy;
		gvy3=gvy2; gvy2=gvy1;
		
		gvy1 = mass(fvy,Bmd,Jr1d,Js1d) - advect(vy,vx,vy,Bmd,Irm1,Ism1...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bvy =       a(1)*gvy1+a(2)*gvy2+a(3)*gvy3;
		bvy = bvy - mass((b(2)*vy1+b(3)*vy2+b(4)*vy3),Bmd,Jr1d,Js1d);
		bvy = bvy - hmhltz(vyb,visc0,b(1),Bmd,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
		bvy = bvy + py;
		bvy = ABu(Ryvy,Rxvy,bvy);

		vyh = visc_slv(bvy,Srvy,Ssvy,Livy,slv);
		vy  = ABu(Ryvy',Rxvy',vyh) + vyb;

		% pressure projection
		if(ifpres)
 			[vx,vy,pr] = pres_proj(vx,vy,pr1...
						,b(1),Bim1,Rxvx,Ryvx,Rxvy,Ryvy,slv...
						,Bm2,Jr12,Js12,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1...
						,Srpr,Sspr,Lipr);
		end
		
	end
	if(ifps)
		 ps3= ps2;  ps2= ps1;  ps1 = ps;
		gps3=gps2; gps2=gps1;
		
		gps1 = mass(fps,Bmd,Jr1d,Js1d) - advect(ps,vx,vy,Bmd,Irm1,Ism1...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bps =       a(1)*gps1+a(2)*gps2+a(3)*gps3;
		bps = bps - mass((b(2)*ps1+b(3)*ps2+b(4)*ps3),Bmd,Jr1d,Js1d);
		bps = bps - hmhltz(psb,visc1,b(1),Bmd,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
		bps = ABu(Ryps,Rxps,bps);

		psh = visc_slv(bps,Srps,Ssps,Lips,slv);
		ps  = ABu(Ryps',Rxps',psh) + psb;

	end

	%[dot(Bm1,vx),dot(Bm1,vy),dot(Bm2,pr),dot(Bm1,ps)]
	if(mod(it,50)==0)

		% log
		[dot(Bm1,vx),dot(Bm1,vy),dot(Bm2,pr),dot(Bm1,ps)]

		if(ifkov)
			[dot(vxe-vx,Bm1),dot(vye-vy,Bm1)]
		end

		% vis
		%surf(xm1,ym1,ps); grid on;

		%omega = vort(vx,vy,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
		%surf(xm1,ym1,omega); grid on;

	  	%quiver(xm1,ym1,vx,vy);
		%contour(xm1,ym1,vx,100);
		surf(xm1,ym1,vx); shading interp;

	   	title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);
		view(2)

		pause(0.01)
	end

	if(blowup(vx,vy,pr,ps));it, return; end;

end
%-------------------------------------------------------------------------------
% post process

['Finished Timestepping']

%surf(xm1,ym1,ps-pse); grid on;
surf(xm1,ym1,vx); shading interp;
%surf(xm1,ym1,ps); grid on;
%quiver(xm1,ym1,vx,vy); grid on;
title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);
view(2)

[dot(vx,Bm1.*vx),dot(vy,Bm1.*vy),dot(pr,Bm2.*pr),dot(ps,Bm1.*ps)]

%===============================================================================
end % driver
%===============================================================================
