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
%	- periodic BC
%
%-------------------------------------------------------------------------------

clf; format compact; format shorte;

n1 = 32;
nd = ceil(1.5*n1);

[z1,w1] = zwgll(n1-1);
[zd,wd] = zwgll(nd-1);

I1 = eye(n1);
Id = eye(nd);

%-------------------------------------------------------------------------------
% geometry

Lx = 1;

[x1] = Lx*ndgrid(z1);
[xd] = Lx*ndgrid(zd);

%-------------------------------------------------------------------------------
% data

% diffusivity
visc = 1e-3;

% initial condition
u = 0*x1; u(end)=1;

% velocity
v = 0*x1;

% forcing
f = 0*x1;
f = sin(pi*x1);

% BC
ub = u;

% Restrictions
R = I1(2:end-1,:);   % dir-dir

ifperiodic = 0;

% T=0 ==> steady
T   = 0.0;
CFL = 0.5;

%------------------------------------------------------------------------------
% setup

% time stepper
dx = min(diff(x1));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0); nt=1;dt=0; end; % steady

% diff matrix
Dr1 = dhat(z1);
Drd = dhat(zd);

D1 = (1/Lx)*Dr1;
Dd = (1/Lx)*Drd;

% interp matrix
J1d = interp_mat(zd,z1); % n1 to nd

% mass matrices
B1  = Lx*diag(w1);
Bd  = Lx*diag(wd);
B1i = (1/Lx)*diag(1./w1);

% mask
msk = diag(R'*R);

% solve
A1 = D1'*B1*D1;

B = R*B1*R';
A = R*A1*R';

%------------------------------------------------------------------------------
% time advance

% initialize time
time = 0;

% initialize histories
time0 = 0;
time1 = 0;
time2 = 0;

u0 = 0*u;
u1 = u0;
u2 = u0;
g1 = 0*R*u;
g2 = g1;

for it=1:nt

	time3=time2; time2=time1; time1 = time;
	time = time + dt;

	if(it<=3)
		[a,b] = bdfext3([time time1 time2 time3]);
		if(T==0) a=0*a; b=0*b; a(1)=1; end; % steady
		H = R*(b(1)*B1 + visc*A1)*R';
	end;

	u3=u2; u2=u1; u1= u;
	g3=g2; g2=g1;
	
    g1 = R*mass(f,B1) - R*advect(u,v,Bd,D1,J1d);

	r =      a(1)*g1+a(2)*g2+a(3)*g3;
	r = r - R*mass(b(2)*u1+b(3)*u2+b(4)*u3,B1);
	r = r - R*hmhltz(ub,b(1),visc,B1,D1);

	uh = H \ r;

	u  = R'*uh + ub;

	% vis
	if(mod(it,100)==0)
		plot(x1,u,'linewidth',2.0);
	   	title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);
		pause(0.01)
	end

	if(blowup(u)) return; end;

end
%-------------------------------------------------------------------------------
% post process

['Finished Timestepping']

plot(x1,u,'linewidth',2.0);
title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);

%===============================================================================
end % driver
%===============================================================================
