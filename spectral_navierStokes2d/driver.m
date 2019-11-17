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
%	- Lid driven cavity failing as nx1/ny1 > 40
%	- verify possion, heat equation, and advection diffusion
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
visc1 = 1e-3;

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
Ryps = Ism1(2:end-1,:);                  % dir-dir

ifxperiodic = 0;
ifyperiodic = 0;

ifvel  = 0;    % evolve velocity field
ifps   = 1;    % evolve passive scalar per advection diffusion
ifconv = 0;    % advect fields per (vx,vy)
ifpres = 0;    % project velocity field onto a div-free subspace

% T=0 ==> steady
T   = 1.0;
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
Bm1 = Jm1.*(wrm1*wsm1');
Bm2 = Jm2.*(wrm2*wsm2');
Bmd = Jmd.*(wrmd*wsmd');
Bm1i= 1 ./ Bm1;
Bm2i= 1 ./ Bm2;

% mask
mskvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)';
mskvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
mskps = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';

% laplace operator setup
g11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
g12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
g22 = Bmd .* (sxmd.*sxmd + symd.*symd);

if(slv==1) % fast diagonalization setup

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
	
	[Srvx,Lrvx] = eig(Arvx,Brvx); Srivx = inv(Srvx);
	[Ssvx,Lsvx] = eig(Asvx,Bsvx); Ssivx = inv(Ssvx);
	Lvx = visc0 * (diag(Lrvx) + diag(Lsvx)');
	
	[Srvy,Lrvy] = eig(Arvy,Brvy); Srivy = inv(Srvy);
	[Ssvy,Lsvy] = eig(Asvy,Bsvy); Ssivy = inv(Ssvy);
	Lvy = visc0 * (diag(Lrvy) + diag(Lsvy)');

	% Passive Scalar
	
	Brps = Rxps*Brv*Rxps';
	Bsps = Ryps*Bsv*Ryps';
	Arps = Rxps*Arv*Rxps';
	Asps = Ryps*Asv*Ryps';
	
	[Srps,Lrps] = eig(Arps,Brps); Srips = inv(Srps);
	[Ssps,Lsps] = eig(Asps,Bsps); Ssips = inv(Ssps);
	Lps = visc1 * (diag(Lrps) + diag(Lsps)');
	
	% Pressure
	Brpr = (Lx/2)*diag(wrm2);
	Bspr = (Ly/2)*diag(wsm2);
	Brvi = (2/Lx)*diag(1./wrm1);
	Bsvi = (2/Ly)*diag(1./wsm1);
	
	Brpr = Brpr*Jr12*(Drv*Brvi*Drv')*Jr12'*Brpr;
	Arpr = Brpr*Jr12*(    Brvi     )*Jr12'*Brpr;
	
	Bspr = Bspr*Js12*(Dsv*Bsvi*Dsv')*Js12'*Bspr;
	Aspr = Bspr*Js12*(    Bsvi     )*Js12'*Bspr;

	[Srpr,Lrpr] = eig(Arpr,Brpr); Sripr = inv(Srpr);
	[Sspr,Lspr] = eig(Aspr,Bspr); Ssipr = inv(Sspr);
	Lpr = diag(Lrpr) + diag(Lspr)';

end

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
fvx1 = vx0;
fvx2 = vx0;
% vy
vy0  = vx0;
vy1  = vx0;
vy2  = vx0;
fvy1 = vx0;
fvy2 = vx0;
% ps
ps0  = vx0;
ps1  = vx0;
ps2  = vx0;
fps1 = vx0;
fps2 = vx0;

for it=1:nt

	time3=time2; time2=time1; time1 = time;
	time = time + dt;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(T  ==0) a=0*a; b=0*b; a(1)=1;   end; % steady
		if(slv==1) Lvxi = 1    ./ (b(1) + Lvx); % FDM
		           Lvyi = 1    ./ (b(1) + Lvy);
		           Lpsi = 1    ./ (b(1) + Lps);
		           Lpri = b(1) ./ (       Lpr); end;
	end;

	if(ifvel)
		vx3=vx2; vx2=vx1; vx1 = vx;
		vy3=vy2; vy2=vy1; vy1 = vy;
						  pr1 = pr;

		% pressure forcing
		if(ifpres) [pr1x,pr1y] = vgradp(pr1,mskvx,mskvy,Bm2,Jr12,Js12,Irm1...
									   ,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
		else pr1x = 0*xm1; pr1y = 0*xm1;
		end
		
		fvx3=fvx2; fvx2=fvx1; fvx1=bdf_expl(vx1,vxb,visc0,mskvx,fvx,vx,vy)+pr1x;
		fvy3=fvy2; fvy2=fvy1; fvy1=bdf_expl(vy1,vyb,visc0,mskvy,fvy,vx,vy)+pr1y;

		bvx = a(1)*fvx1+a(2)*fvx2+a(3)*fvx3 - Bm1.*(b(2)*vx1+b(3)*vx2+b(4)*vx3);
		bvy = a(1)*fvy1+a(2)*fvy2+a(3)*fvy3 - Bm1.*(b(2)*vy1+b(3)*vy2+b(4)*vy3);

		vxh = hmhltz_slv(bvx,mskvx,Bm1i,Srvx,Ssvx,Srivx,Ssivx,Rxvx,Ryvx,Lvxi,slv);
		vyh = hmhltz_slv(bvy,mskvy,Bm1i,Srvy,Ssvy,Srivy,Ssivy,Rxvy,Ryvy,Lvyi,slv);

		vx  = vxh + vxb;
		vy  = vyh + vyb;

		if(ifpres); [vx,vy,pr] = pres_proj(vx,vy,pr1); end;
	end
	if(ifps)
		 ps3= ps2;  ps2= ps1;  ps1 = ps;
		fps3=fps2; fps2=fps1; fps1 = bdf_expl(ps1,psb,visc1,mskps,fps,vx,vy);

		bps = a(1)*fps1+a(2)*fps2+a(3)*fps3 - Bm1.*(b(2)*ps1+b(3)*ps2+b(4)*ps3);

		psh = hmhltz_slv(bps,mskps,Bm1i,Srps,Ssps,Srips,Ssips,Rxps,Ryps,Lpsi,slv);
		ps  = psh + psb;
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
%
%	Helper functions
%
%===============================================================================
% BDF - implicit OP
function [Hu] =  hmhltz(uhm,mskhm,vischm)
	Hu =    vischm*lapl(uhm,mskhm,Jr1d,Js1d,Drm1,Dsm1,g11,g12,g22);
	Hu = Hu + b(1)*mass(uhm,mskhm,Bmd,Jr1d,Js1d);
end

%-------------------------------------------------------------------------------
% BDF - explicit OP
function [Fu] = bdf_expl(uexp,ubexp,viscexp,mskexp,fexp,cx,cy)
	Fu = mass(fexp,mskexp,Bmd,Jr1d,Js1d);                       % forcing
	if(ifconv);
		Fu = Fu -advect(uexp,mskexp,cx,cy,Bmd,Irm1,Ism1,Jr1d... % convection
				       ,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);
	end
	Fu = Fu - hmhltz(ubexp,1+0*mskexp,viscexp);                 % dirichlet BC
end

%-------------------------------------------------------------------------------
% pressure project

function [ux,uy,p] = pres_proj(cx,cy,pr_prev)

	g = -qdivu(cx,cy,mskvx,mskvy,Bm2,Jr12,Js12,Irm1...
			  ,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

	if(slv==1) % FDM
		p = fdm(g,Bm2i,Srpr,Sspr,Sripr,Ssipr,Irm2,Ism2,Lpri);
	end

	[px,py] = vgradp(p,mskvx,mskvy,Bm2,Jr12,Js12...
	                ,Irm1,Ism1,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

	ux = cx + 1/(b(1))*Bm1i.*px;
	uy = cy + 1/(b(1))*Bm1i.*py;

	p = p + pr_prev;

end
%===============================================================================
end % driver
%===============================================================================
