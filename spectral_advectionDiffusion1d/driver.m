%===============================================================================
%
%	Driver function for 1 D diffusion equation
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
%
%-------------------------------------------------------------------------------

clf; format compact; format shorte;

n1 = 32;
nd = ceil(1.5*n1);

[z1,w1] = zwgll(n1-1);
[zd,wd] = zwgll(nd-1);

D1 = dhat(z1);
Dd = dhat(zd);

I1 = eye(n1);
Id = eye(nd);

J1d = interp_mat(zd,z1); % n1 to nd

%-------------------------------------------------------------------------------
% geometry

[x1] = ndgrid(z1);
[xd] = ndgrid(zd);

%-------------------------------------------------------------------------------
% data

% solver --> 0: CG, 1: FDM
slv=1;

% viscosity (velocity, passive scalar)
visc0 = 1/2e1;
visc1 = 1e-3;

% initial condition
v  = 0*xm1;
s  = 0*xm1; s(:,end)=1;

% forcing
fv = 0*x1;
fs = 0*x1;

% BC
vb = v;
sb = s;

% Restrictions
Rv = I1(2:end-1,:); % vx             % dir-dir
Rs = I1(2:end-1,:); % ps             % dir-dir

ifxperiodic = 0;

ifvel  = 0;    % evolve velocity field
ifs    = 1;    % evolve passive scalar per advection diffusion
ifconv = 0;    % advect fields per (vx,vy)

% T=0 ==> steady
T   = 1.0;
CFL = 0.5;

%------------------------------------------------------------------------------
% setup

% time stepper
dx = min(min(diff(x1)));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0); nt=1;dt=0; end; % steady

% mass
B1 = w1;
Bd = wd;
B1i= 1 ./ B1;

% mask
mskv = diag(Rv'*Rv);
msks = diag(Rs'*Rs);

% direct solve

A1 = D1'*B1*D1;

Bv = Rv*B1*Rv';
Av = Rv*A1*Rv';

Bs = Rs*B1*Rs';
As = Rs*A1*Rs';

[Sv,Lv] = eig(Av,Bv); Siv = inv(Sv);
Lv = visc0 * diag(Lv);

[Ss,Ls] = eig(As,Bs); Sis = inv(Ss);
Ls = visc1 * diag(Ls);

%------------------------------------------------------------------------------
% time advance

time = 0;

% initialize histories
time0 = 0;
time1 = 0;
time2 = 0;
% v
v  = v*0;
v  = v;
v  = v;
fv1 = v0;
fv2 = v0;
% s
s0  = v0;
s1  = v0;
s2  = v0;
fs1 = v0;
fs2 = v0;

for it=1:nt

	time3=time2; time2=time1; time1 = time;
	time = time + dt;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(T  ==0) a=0*a; b=0*b; a(1)=1; end; % steady
		if(slv==1) Lvi = 1    ./ (b(1) + Lv);
		           Lsi = 1    ./ (b(1) + Ls); end
	end;

	if(ifvel)
		 v3= v2;  v2= v1;  v1= vx;
		fv3=fv2; fv2=fv1; fv1=bdf_expl(vx1,vxb,visc0,mskvx,fvx,vx,vy);

		bv = a(1)*fv1+a(2)*fv2+a(3)*fv3 - B1.*(b(2)*v1+b(3)*v2+b(4)*v3);
		vh = hmhltz_slv(bv,mskv,Bi,Srvx,Ssvx,Srivx,Ssivx,Rxvx,Ryvx,Lvxi,slv);
		v  = vh + vb;

		if(ifpres); [vx,vy,pr] = pres_proj(vx,vy,pr1); end;
	end
	if(ifps)
		 s3= s2;  s2= s1;  s1 = s;
		fs3=fs2; fs2=fs1; fs1 = bdf_expl(s1,sb,visc1,mskps,fps,v);

		bs = a(1)*fps1+a(2)*fps2+a(3)*fps3 - Bm1.*(b(2)*ps1+b(3)*ps2+b(4)*ps3);

		sh = hmhltz_slv(bs,msks,B1i,Srps,Ssps,Srips,Ssips,Rxps,Ryps,Lpsi,slv);
		s  = sh + sb;
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

	if(blowup(v,s)) return; end;

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

%===============================================================================
end % driver
%===============================================================================
