%===============================================================================
%
%	Driver function for 1 D diffusion equation
%
%	du/dt + \vect{c}\dot\grad{u}  = f + nu*\del^2 u
%
%   + Dirichlet/Neumann BC
%
%===============================================================================
%function driver
%
%-------------------------------------------------------------------------------
%
%	/todo
%
%-------------------------------------------------------------------------------

clf; format compact; format shorte;

E  = 4;
n1 = 8;
nd = ceil(1.5*n1);

[z1,w1] = zwgll(n1-1);
[zd,wd] = zwgll(nd-1);

Q1 = semq(E,n1,0);
Qd = semq(E,nd,0);

n1g = size(Q1,2);
ndg = size(Qd,2);

x1 = semmesh(E,n1,0);
xd = semmesh(E,nd,0);

Dr1 = dhat(z1);
Drd = dhat(zd);

D1 = kron(speye(E),Dr1);
Dd = kron(speye(E),Drd);

Jr1d = interp_mat(zd,z1);
J1d  = kron(speye(E),Jr1d);

I1 = speye(n1g);
Id = speye(ndg);

R = I1(2:end-1,:);

B1  = (2/E) * kron(ones(E,1),w1);
Bd  = (2/E) * kron(ones(E,1),wd);
Bi1 = 1 ./ B1;

%-------------------------------------------------------------------------------
% setup

% mask
msk = diag(R'*R);

if(ifperiodic)
	R = [eye(n1-1),[1;zeros(n1-2,1)]];
end;

% time stepper
dx = min(diff(x1));
dt = dx*CFL/1;
nt = floor(T/dt);
dt = T/nt;

if(T==0); nt=1;dt=0; end; % steady

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
	if(mod(it,10)==0 | time>=T-1e-6)
		up = J1p*u;
		plot(xp,up,'linewidth',2.0);
	   	title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);
		drawnow;
	end

	if(blowup(u)) return; end;

end
%-------------------------------------------------------------------------------
% post process

['Finished Timestepping']

title(['t=',num2str(time),', Step ',num2str(it),' CFL=',num2str(CFL)]);

%===============================================================================
%end % driver
%===============================================================================
