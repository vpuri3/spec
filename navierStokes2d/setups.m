%------------------------------------------------------------------------------
% kovazny flow

a = -0.5; lx = 2.5; ly=2.0;
xx = a + lx/2 * (zrm1+1) ; yy = a + ly/2 * (zsm1+1); [xm1,ym1] = ndgrid(xx,yy);
xx = a + lx/2 * (zrm2+1) ; yy = a + ly/2 * (zsm2+1); [xm2,ym2] = ndgrid(xx,yy);
xx = a + lx/2 * (zrmd+1) ; yy = a + ly/2 * (zsmd+1); [xmd,ymd] = ndgrid(xx,yy);  

Re = 40;  % target Re
Re = 10;
nu = 1/Re;

% boundary conditions
vx = 0*xm1;
vy = 0*xm1;

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

pr = 0*xm2; % no pressure BCs in this case
fx = 0*xm1; % forcing
fy = 0*xm1;

% BC --> smooth functinos
vxb = vx;
vyb = vy;
end

ifvel  = 1;
ifconv = 1;
ifpres = 0;

%-------------------------------------------------------------------------------
% lid driven cavity

% solver --> 0: CG, 1: FDM
slv=1;

% viscosity (velocity, passive scalar)
visc0 = 1e-2;
visc1 = 1e-0;

% initial condition
vx  = 0*xm1; vx(:,end)=1;
vy  = 0*xm1;
ps  = 0*xm1;
pr  = 0*xm2;

% forcing
fvx = 0*xm1;
fvy = 0*xm1;
fps = 0*xm1; fps = sin(pi*xm1).*sin(pi*ym1); pse = fps/2/pi/pi/visc1;

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
T   = 2e1;
CFL = 0.5;

%------------------------------------------------------------------------------

