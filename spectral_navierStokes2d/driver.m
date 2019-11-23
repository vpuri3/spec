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
%	- periodicity
%	- testing
%
%-------------------------------------------------------------------------------

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

%-------------------------------------------------------------------------------
% geometry

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xm2,ym2] = ndgrid(zrm2,zsm2);
[xmd,ymd] = ndgrid(zrmd,zsmd);

%-------------------------------------------------------------------------------
% data

% solver --> 0: CG, 1: FDM
slv=1;

% viscosity (velocity, passive scalar)
visc0 = 1/2e1;
visc1 = 1e-2;

% initial condition
vx  = 0*xm1;
vy  = 0*xm1;
ps  = 0*xm1; ps(:,end)=1;
pr  = 0*xm2;

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
Rxps = Irm1(1:end-0,:); % ps             % neu-neu
Ryps = Ism1(2:end-1,:);                  % dir-dir

ifxperiodic = 0;
ifyperiodic = 0;

ifvel  = 0;    % evolve velocity field
ifps   = 1;    % evolve passive scalar per advection diffusion
ifpres = 1;    % project velocity field onto a div-free subspace

% T=0 ==> steady
T   = 5.0;
CFL = 0.5;

%------------------------------------------------------------------------------
% setup

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
Bm1i = 1   ./ Bm1;
Bm2i = 1   ./ Bm2;

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
mskps = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';

% laplace operator setup
g11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
g12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
g22 = Bmd .* (sxmd.*sxmd + symd.*symd);

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

RBivx = 1 ./ ABu(Ryvx,Rxvx,Bm1);
[Srvx,Lrvx] = eig(Arvx,Brvx); Srivx = inv(Srvx);
[Ssvx,Lsvx] = eig(Asvx,Bsvx); Ssivx = inv(Ssvx);
Lvx = visc0 * (diag(Lrvx) + diag(Lsvx)');

RBivy = 1 ./ ABu(Ryvy,Rxvy,Bm1);
[Srvy,Lrvy] = eig(Arvy,Brvy); Srivy = inv(Srvy);
[Ssvy,Lsvy] = eig(Asvy,Bsvy); Ssivy = inv(Ssvy);
Lvy = visc0 * (diag(Lrvy) + diag(Lsvy)');

% Passive Scalar

Brps = Rxps*Brv*Rxps';
Bsps = Ryps*Bsv*Ryps';
Arps = Rxps*Arv*Rxps';
Asps = Ryps*Asv*Ryps';

RBips = 1 ./ ABu(Ryps,Rxps,Bm1);
[Srps,Lrps] = eig(Arps,Brps); Srips = inv(Srps);
[Ssps,Lsps] = eig(Asps,Bsps); Ssips = inv(Ssps);
Lps = visc1 * (diag(Lrps) + diag(Lsps)');

% Pressure
Brpr = (Lx/2)*diag(wrm2);
Bspr = (Ly/2)*diag(wsm2);
Briv = (2/Lx)*diag(1./wrm1);
Bsiv = (2/Ly)*diag(1./wsm1);

Brpr = Brpr*Jr12*(Drv*Briv*Drv')*Jr12'*Brpr;
Arpr = Brpr*Jr12*(    Briv     )*Jr12'*Brpr;

Bspr = Bspr*Js12*(Dsv*Bsiv*Dsv')*Js12'*Bspr;
Aspr = Bspr*Js12*(    Bsiv     )*Js12'*Bspr;

[Srpr,Lrpr] = eig(Arpr,Brpr); Sripr = inv(Srpr);
[Sspr,Lspr] = eig(Aspr,Bspr); Ssipr = inv(Sspr);
Lpr = diag(Lrpr) + diag(Lspr)';

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
	end;

	if(ifvel)
		% pr 
		if(ifpres)
			pr1 = pr;
			[px,py] = vgradp(pr1,Bm2,Jr12,Js12,Irm1,Ism1...
	 						,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
		else
			px = 0*vx;
			py = 0*vy;
		end

		% vx
		 vx3= vx2;  vx2= vx1;  vx1 = vx;
		gvx3=gvx2; gvx2=gvx1;
		
		gvx1 = mass(fvx,Bmd,Jr1d,Js1d) - advect(vx,vx,vy,Bmd,Irm1,Ism1...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bvx =       a(1)*gvx1+a(2)*gvx2+a(3)*gvx3;
		bvx = bvx - mass((b(2)*vx1+b(3)*vx2+b(4)*vx3),Bmd,Jr1d,Js1d);
		bvx = bvx - hmhltz(vxb,visc1,b(1),Bmd,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
		bvx = bvx + px;
		bvx = ABu(Ryvx,Rxvx,bvx);

		vxh = visc_slv(bvx,RBivx,Srvx,Ssvx,Srivx,Ssivx,Livx,slv);

		vx  = ABu(Ryvx',Rxvx',vxh) + vxb;

		% vy
		 vy3= vy2;  vy2= vy1;  vy1 = vy;
		gvy3=gvy2; gvy2=gvy1;
		
		gvy1 = mass(fvy,Bmd,Jr1d,Js1d) - advect(vy,vx,vy,Bmd,Irm1,Ism1...
									  ,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

		bvy =       a(1)*gvy1+a(2)*gvy2+a(3)*gvy3;
		bvy = bvy - mass((b(2)*vy1+b(3)*vy2+b(4)*vy3),Bmd,Jr1d,Js1d);
		bvy = bvy - hmhltz(vyb,visc1,b(1),Bmd,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
		bvx = bvx + py;
		bvy = ABu(Ryvy,Rxvy,bvy);

		vyh = visc_slv(bvy,RBivy,Srvy,Ssvy,Srivy,Ssivy,Livy,slv);

		vy  = ABu(Ryvy',Rxvy',vyh) + vyb;

		% pressure
		if(ifpres)
 			[vx,vy,pr] = pres_proj(vx,vy,pr1...
						,Bm2,Jr12,Js12,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1...
						,Bm2i,Srpr,Sspr,Sripr,Ssipr,Lipr)
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

		psh = visc_slv(bps,RBips,Srps,Ssps,Srips,Ssips,Lips,slv);

		ps  = ABu(Ryps',Rxps',psh) + psb;

	end

	% vis
	if(mod(it,100)==0)
		% pseudocolor subplots for viewing velocity field
		%omega = vort(vx,vy,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	  	%quiver(xm1,ym1,vx,vy); grid on;
		%[dot(vx,Bm1.*vx),dot(vy,Bm1.*vy)]
		surf(xm1,ym1,ps); shading interp
	   	title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);
		pause(0.01)
	end

	if(blowup(vx,vy,pr,ps)) return; end;

end
%-------------------------------------------------------------------------------
% post process

['Finished Timestepping']

surf(xm1,ym1,ps); shading interp
title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);

%===============================================================================
end % driver
%===============================================================================
