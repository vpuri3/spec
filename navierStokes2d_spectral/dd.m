%------------
ifannulus=1;
ifwavy   =0;
ifrb     =1;

nx1 = 32;
ny1 = 32;

slv=0; % 0: CG, 1: FDM

%-------------------------------------------------------------------------------
% geometry

%[xm1,ym1] = ndgrid(zrm1,zsm1);
%[xm2,ym2] = ndgrid(zrm2,zsm2);
%[xmd,ymd] = ndgrid(zrmd,zsmd);
%[xmp,ymp] = ndgrid(zrmp,zsmp);
%
kc = 3.117; a = 0; Lx = 2*pi/kc; Ly=1;
xx = a + Lx/2 * (zrm1+1); yy = a + Ly/2 * (zsm1+1); [xm1,ym1]=ndgrid(xx,yy);
xx = a + Lx/2 * (zrm2+1); yy = a + Ly/2 * (zsm2+1); [xm2,ym2]=ndgrid(xx,yy);
xx = a + Lx/2 * (zrmd+1); yy = a + Ly/2 * (zsmd+1); [xmd,ymd]=ndgrid(xx,yy);
xx = a + Lx/2 * (zrmp+1); yy = a + Ly/2 * (zsmp+1); [xmp,ymp]=ndgrid(xx,yy);

if(ifwavy)
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = wwall(zrm1,zsm1);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = wwall(zrm2,zsm2);
[xm2,ym2] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm2,zsm2);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = wwall(zrmd,zsmd);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = wwall(zrmp,zsmp);
[xmp,ymp] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmp,zsmp);
end

if(ifannulus)
r0=1; r1=2; span=2*pi;
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrm1,zsm1);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrm2,zsm2);
[xm2,ym2] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm2,zsm2);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrmd,zsmd);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zrmp,zsmp);
[xmp,ymp] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmp,zsmp);
end

if(ifannulus)
	Lx = 0.5*span*(r1*r1-r0*r0);
	Ly = r1-r0;
else
	Lx = abs(sum(xm1(end,:)-xm1(1,:))) / nx1;
	Ly = abs(sum(ym1(:,end)-ym1(:,1))) / ny1;
end

%------------------------------------------------------------------------------
casename = 'RB'; cname = 'rb';

Ra = 100; Pr=1;

visc0 = Pr;
visc1 = 1;

vx  = xm1*0;
vy  = ym1*0;
pr  = 0*xm2;
ps  = 0*xm1; ps(:,1)=1;

fvx = 0*xm1; fvy = 0*xm1; fps = 0*xm1; 
vxb = vx; vyb = vy; psb = ps;

Rxvx = Irm1(2:end-1,:); Ryvx = Irm1(2:end-1,:);
Rxvy = Ism1(2:end-1,:); Ryvy = Ism1(2:end-1,:);
Rxps = Irm1(2:end-1,:); Ryps = Ism1(2:end-1,:);

ifxperiodic = 0;
ifyperiodic = 0;

T   = 50.0;
CFL = 0.5;

%------------------------------------------------------------------------------
