%
% kovazny

casename = 'Kovasznay Flow'; cname = 'kov';

% viscosity (velocity, passive scalar)
Re = 40;
visc0 = 1/Re;
visc1 = 1e-0;

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
T   = 10.0;
CFL = 0.1;

%--- diffusion check
%visc1 = 1e-0;
%fps=1+0*xm1g; T=10;

%--- convection check
%T = 2*pi;
%CFL = 0.3;
%visc1 = 0e-0;
%vx= ym1g;
%vy=-xm1g;
%d2=(xm1g+0.3).^2 + (ym1g-0.0).^2;
%ps=exp(-d2/0.016);
%psb = ps;
