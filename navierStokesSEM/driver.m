%===============================================================================
%
%	Driver function for Navier Stokes solver
%
%	\partial_t u + (u\dot\grad) u = -\grad p + nu*\del^2 u + f
%   				  \grad\dot u = 0
%
%   + Dirichlet/Neumann/Periodic BC
%
%===============================================================================
%
%function driver
%
%-------------------------------------------------------------------------------
%
%	/todo
%	- Preconditioners
%		- Hlmhltz -> add viscous information
%		- Pres    -> Schwarz?
%	- artificial viscosity, stabalization
%
%-------------------------------------------------------------------------------
clear;
clf; fig=gcf;
format compact; format shorte;

% initialization
[ifvel,ifpres,ifps,ifadv,Ex,Ey,nx1,ny1,nx2,ny2,nxd,nyd] = usrinit;

% element operators
[zrm1,wrm1] = zwgll(nx1-1); [zsm1,wsm1] = zwgll(ny1-1); % vel, scalar
[zrm2,wrm2] = zwgll(nx2-1); [zsm2,wsm2] = zwgll(ny2-1); % pres
[zrmd,wrmd] = zwgll(nxd-1); [zsmd,wsmd] = zwgll(nyd-1); % dealias

Drm1 = dhat(zrm1); Dsm1 = dhat(zsm1);
Drm2 = dhat(zrm2); Dsm2 = dhat(zsm2);
Drmd = dhat(zrmd); Dsmd = dhat(zsmd);

Jr1d = interp_mat(zrmd,zrm1); Js1d = interp_mat(zsmd,zsm1); % vel  -> dealias
Jr21 = interp_mat(zrm1,zrm2); Js21 = interp_mat(zsm1,zsm2); % pres -> vel

% local operators
wxm1 = kron(ones(Ex,1),wrm1); wym1 = kron(ones(Ey,1),wsm1);
wxm2 = kron(ones(Ex,1),wrm2); wym2 = kron(ones(Ey,1),wsm2);
wxmd = kron(ones(Ex,1),wrmd); wymd = kron(ones(Ey,1),wsmd);

Dxm1 = kron(speye(Ex),Drm1); Dym1 = kron(speye(Ey),Dsm1);
Dxm2 = kron(speye(Ex),Drm2); Dym2 = kron(speye(Ey),Dsm2);
Dxmd = kron(speye(Ex),Drmd); Dymd = kron(speye(Ey),Dsmd);

Jx1d = kron(speye(Ex),Jr1d); Jy1d = kron(speye(Ey),Js1d);
Jx21 = kron(speye(Ex),Jr21); Jy21 = kron(speye(Ey),Js21);

% bc
[ifxprdc,ifyprdc,Mvx,Mvy,Mps] = usrbc(Ex,Ey,nx1,ny1);

% global -> local operator
Qx1 = semq(Ex,nx1-1,ifxprdc); Qy1 = semq(Ey,ny1-1,ifyprdc);
Qx2 = semq(Ex,nx2-1,ifxprdc); Qy2 = semq(Ey,ny2-1,ifyprdc);

% geom
[xm1,ym1] = usrgeom(nx1,ny1,Ex,Ey);
[xm2,ym2] = usrgeom(nx2,ny2,Ex,Ey);
[xmd,ymd] = usrgeom(nxd,nyd,Ex,Ey);

% case setup, ic, dt
[visc0,visc1,T,time,vx,vy,pr,ps,dt,nstep] = usrcase(xm1,ym1,xm2,ym2);

% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Dxm1,Dym1);
[Jm2,Jim2,rxm2,rym2,sxm2,sym2] = jac2d(xm2,ym2,Dxm2,Dym2);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Dxmd,Dymd);

% diag mass mats
Bm1  = Jm1 .* (wxm1*wym1');
Bm2  = Jm2 .* (wxm2*wym2');
Bmd  = Jmd .* (wxmd*wymd');
Bim1 = 1   ./ Bm1;

vol = dot(Bm1,1+0*Bm1);

% lapl op setup
g11 = Bm1 .* (rxm1 .* rxm1 + rym1 .* rym1);
g12 = Bm1 .* (rxm1 .* sxm1 + rym1 .* sym1);
g22 = Bm1 .* (sxm1 .* sxm1 + sym1 .* sym1);

%------------------------------------------------------------------------------
% time advance

% initialize histories
time1 = time; time2 = time; time3=time;
vx1 = vx*0; vx2 = vx1; vx3 = vx2; gvx1 = vx1; gvx2 = vx1; gvx3 = vx1;
vy1 = vy*0; vy2 = vy1; vy3 = vy2; gvy1 = vy1; gvy2 = vy1; gvy3 = vy1;
ps1 = ps*0; ps2 = ps1; ps3 = ps2; gps1 = ps1; gps2 = ps1; gps3 = ps1;
pr1 = pr*0; pr2 = pr*0;

% vis, log
usrchk(xm1,ym1,xm2,ym2,vx,vy,pr,ps,time,T,0,nstep,fig);

for istep=1:nstep

	% update histories
	time3=time2; time2=time1; time1=time;

	vx3=vx2; vx2=vx1; vx1=vx; gvx3=gvx2; gvx2=gvx1;
	vy3=vy2; vy2=vy1; vy1=vy; gvy3=gvy2; gvy2=gvy1;
	ps3=ps2; ps2=ps1; ps1=ps; gps3=gps2; gps2=gps1;
			 pr2=pr1; pr1=pr;
	
	% update time
	time = time + dt;

	% update BC, forcing
	[vxb,vyb,prb,psb,fvx,fvy,fps] = usrf(xm1,ym1,xm2,ym2,time);

	% update time-stepper
	if(istep<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(istep==0) ap=[1 0]; elseif(istep==1) ap=[2 -1]; end
		if(T==0) a=0*a; b=0*b; a(1)=1;			           end % steady heat eqn
	end

	% solve
	if(ifps)
		gps1 = -ifadv*advect(ps1,vx1,vy1,[],[],[],Bmd...
						    ,Jx1d,Jy1d,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);

		bps =       a(1)*gps1+a(2)*gps2+a(3)*gps3;
		bps = bps + mass(fps,Bm1,[],[],[]);
		bps = bps - mass((b(2)*ps1+b(3)*ps2+b(4)*ps3),Bm1,[],[],[]);
		bps = bps - hlmhltz(psb,visc1,b(1),[],[],[],Bm1,Dxm1,Dym1,g11,g12,g22);
		bps = mass(bps,[],Mps,Qx1,Qy1);

		psh = pcg_visc(bps,visc1,b(1),Mps,Qx1,Qy1...
					  ,Bm1,Bim1,Dxm1,Dym1,g11,g12,g22,1e-8,1e3);

		ps  = mask(psh,Mps) + psb;
	end

	if(ifvel)
		
		gvx1 = -ifadv*advect(vx1,vx1,vy1,[],[],[],Bmd...
			                ,Jx1d,Jy1d,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);
		gvy1 = -ifadv*advect(vy1,vx1,vy1,[],[],[],Bmd...
						    ,Jx1d,Jy1d,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);

		% pressure forcing
		pr = ap(1)*pr1 + ap(2)*pr2;
		if(ifpres) [px,py]=gradp(pr+prb,[],[]...
					  			 ,Bm1,Jx21,Jy21,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);
		else px=0*xm1; py=0*xm1;
		end

		% viscous solve
		bvx =       a(1)*gvx1+a(2)*gvx2+a(3)*gvx3;
		bvx = bvx + mass(fvx,Bm1,[],[],[]);
		bvx = bvx - mass((b(2)*vx1+b(3)*vx2+b(4)*vx3),Bm1,[],[],[]);
		bvx = bvx - hlmhltz(vxb,visc0,b(1),[],[],[],Bm1,Dxm1,Dym1,g11,g12,g22);
		bvx = bvx + px;
		bvx = mass(bvx,[],Mvx,Qx1,Qy1);

		bvy =       a(1)*gvy1+a(2)*gvy2+a(3)*gvy3;
		bvy = bvy + mass(fvy,Bm1,[],[],[]);
		bvy = bvy - mass((b(2)*vy1+b(3)*vy2+b(4)*vy3),Bm1,[],[],[]);
		bvy = bvy - hlmhltz(vyb,visc0,b(1),[],[],[],Bm1,Dxm1,Dym1,g11,g12,g22);
		bvy = bvy + py;
		bvy = mass(bvy,[],Mvy,Qx1,Qy1);

		vxh = pcg_visc(bvx,visc0,b(1),Mvx,Qx1,Qy1...
					  ,Bm1,Bim1,Dxm1,Dym1,g11,g12,g22,1e-8,1e3);
		vyh = pcg_visc(bvy,visc0,b(1),Mvy,Qx1,Qy1...
					  ,Bm1,Bim1,Dxm1,Dym1,g11,g12,g22,1e-8,1e3);

		vx  = mask(vxh,Mvx) + vxb;
		vy  = mask(vyh,Mvy) + vyb;

		% pressure projection
		if(ifpres)
 			[vx,vy,pr] = pres_proj(vx,vy,pr,b(1),Mvx,Mvy,Qx1,Qy1,Qx2,Qy2...
						,Bm1,Bim1,Jx21,Jy21,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);
		end

	end

	% vis, log
	usrchk(xm1,ym1,xm2,ym2,vx,vy,pr,ps,time,T,istep,nstep,fig);

end

%===============================================================================
%
%	case setup
%
%===============================================================================
function [ifvel,ifpres,ifps,ifadv,Ex,Ey,nx1,ny1,nx2,ny2,nxd,nyd] = usrinit

ifvel  = 1;    % evolve  vel field per NS eqn
ifadv  = 1;    % advect  vel, sclr
ifpres = 1;    % project vel onto a div-free subspace
ifps   = 0;    % evolve sclr per advection diffusion eqn

Ex  = 4;
Ey  = 4;
nx1 = 8;
ny1 = nx1;
nx2 = nx1-2;
ny2 = ny1-2;
nxd = ceil(1.5*nx1); nxd = nxd + rem(nxd,2);
nyd = ceil(1.5*ny1); nyd = nyd + rem(nyd,2);

end
%------------------------------------------------------------------------------
function [ifxprdc,ifyprdc,Mvx,Mvy,Mps] = usrbc(Ex,Ey,nx1,ny1)

ifxprdc = 0;
ifyprdc = 0;

Ixm1 = speye(Ex*nx1);
Iym1 = speye(Ey*ny1);

Rxvx = Ixm1(2:end-1,:); Ryvx = Iym1(2:end-1,:); % restriction
Rxvy = Ixm1(2:end-1,:); Ryvy = Iym1(2:end-1,:);
Rxps = Ixm1(2:end-1,:); Ryps = Iym1(2:end-1,:);

if(ifxprdc) Rxvx = Ixm1; Rxvy = Ixm1; Rxps = Ixm1; end;
if(ifyprdc) Ryvx = Iym1; Ryvy = Iym1; Ryps = Iym1; end;

Mvx = diag(Rxvx'*Rxvx) * diag(Ryvx'*Ryvx)'; % mask
Mvy = diag(Rxvy'*Rxvy) * diag(Ryvy'*Ryvy)';
Mps = diag(Rxps'*Rxps) * diag(Ryps'*Ryps)';

end
%------------------------------------------------------------------------------
function [x,y] = usrgeom(nx,ny,Ex,Ey)

[x1d,~] = semmesh(Ex,nx,0);
[y1d,~] = semmesh(Ey,ny,0);

a = -0.5; lx = 2.5; ly=2.0;%a = -0.5; lx = 1.0; ly=1.0;
xx=a+lx/2*(x1d+1); yy=a+ly/2*(y1d+1); [x,y] = ndgrid(xx,yy);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(x1d,y1d);
[x,y] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,x1d,y1d);

end
%------------------------------------------------------------------------------
function [visc0,visc1,T,time,vx,vy,pr,ps,dt,nstep] = usrcase(xm1,ym1,xm2,ym2)

% visc (vel, ps)
Re = 40;
visc0 = 1/Re;
visc1 = 1e-0;

% initial condition
time = 0;
vx = 0*xm1;
vy = 0*xm1;
ps = 0*xm1;
pr = 0*xm2;

T   = 10.0;
CFL = 0.1;

% time stepper
%xm1g = unique(xm1); dx = min(min(abs(diff(xm1g))));
%ym1g = unique(ym1); dy = min(min(abs(diff(ym1g))));
%dx = min(dx,dy);
%dt = dx*CFL/1;
dt=1e-3;
nstep = floor(T/dt);
dt = T/nstep;

if(T==0) nstep=1; dt=0; end % steady diffusion equation

end
%------------------------------------------------------------------------------
function [vxb,vyb,prb,psb,fvx,fvy,fps] = usrf(xm1,ym1,xm2,ym2,time)

Re = 40;

[vxe,vye] = kov_ex(xm1,ym1,Re);%mesh(xm1,ym1,vxe);pause;
vxb = vxe;
vyb = vye;
prb = 0*xm2;
psb = 0*xm1;

fvx = 0*xm1;
fvy = 0*xm1;
fps = 0*xm1;

end
%------------------------------------------------------------------------------
function usrchk(xm1,ym1,xm2,ym2,vx,vy,pr,ps,time,T,istep,nstep,fig)

blowup(vx,vy,pr,ps,istep);

persistent casename cname mov;

if(istep==0)
	casename = 'Kovazney';
	cname = 'kov';
	mov = []; % movie
end

% vis, log
if(mod(istep,5e1)==0 & istep~=0 | istep==nstep)

	Re = 40;
	[vxe,vye] = kov_ex(xm1,ym1,Re);

	['infty kovazny normalized v-ve']
	[max(max(abs(vx-vxe))),max(max(abs(vy-vye)))] ./...
	[max(max(abs(vxe))),max(max(abs(vye)))]

	% vis
	%om = vort(vx,vy,Qx1,Qy1,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);

	contour(xm1,ym1,vx,20); view(2);colorbar
   	title([casename,', t=',num2str(time,'%4.2f'),' i=',num2str(istep)]);
	drawnow
	mov = [mov,getframe(fig)];
end

if(istep==nstep)
	['Finished Timestepping']
	%['Energy in vx,vy,pr,ps'],[L2(vx,Bm1),L2(vy,Bm1),L2(pr,Bm2),L2(ps,Bm1)]
	
	%% play movie
	%movie(fig,mov,-2,40);
	
	%% save as gif
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
end

end
%===============================================================================
%end % driver
%===============================================================================
