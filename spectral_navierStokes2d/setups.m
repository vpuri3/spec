%-------------------------------------------------------------------------------
% lid driven cavity

slv=1;                           % solver --> 0: CG, 1: FDM, 2: direct

nu = 1/2e1;
vx = 0*xm1; vx(:,end) = 1;
vy = 0*xm1;
pr = 0*xm2;
fx = 0*xm1; % forcing
fy = 0*xm1;

% BC --> smooth functinos
vxb = vx;
vyb = vy;

Rxvx = Irm1(2:end-1,:);                  % dir-dir
Ryvx = Ism1(2:end-1,:);                  % dir-dir
Rxvy = Irm1(2:end-1,:);                  % dir-dir
Ryvy = Ism1(2:end-1,:);                  % dir-dir

ifvel  = 1;
ifconv = 1;
ifpres = 1;

T   = 1; % T=0 ==> steady
CFL = 0.5;

%------------------------------------------------------------------------------
% kovazny flow
ifkov = 1;
if(ifkov);

a = -0.5; lx = 2.5; ly=2.0;
xx = a + lx/2 * (zrm1+1) ; yy = a + ly/2 * (zsm1+1); %linear mapping
[xm1,ym1] = ndgrid(xx,yy);
xx = a + lx/2 * (zrm2+1) ; yy = a + ly/2 * (zsm2+1);
[xm2,ym2] = ndgrid(xx,yy);
xx = a + lx/2 * (zrmd+1) ; yy = a + ly/2 * (zsmd+1);
[xmd,ymd] = ndgrid(xx,yy);  

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

%------------------------------------------------------------------------------

