%
% deform geometry
%
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xm1g,ym1g);
%[xm1g,ym1g] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xm1g,ym1g);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xm2g,ym2g);
%[xm2g,ym2g] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xm2g,ym2g);
%[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(xmdg,ymdg);
%[xmdg,ymdg] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,xmdg,ymdg);

a = -0.5; lx = 2.5; ly=2.0;
a = -1.0; lx = 2.0; ly=2.0;
xx=a+lx/2*(xm1g+1); yy=a+ly/2*(ym1g+1); [xm1g,ym1g] = ndgrid(xx,yy);
xx=a+lx/2*(xm2g+1); yy=a+ly/2*(ym2g+1); [xm2g,ym2g] = ndgrid(xx,yy);
xx=a+lx/2*(xmdg+1); yy=a+ly/2*(ymdg+1); [xmdg,ymdg] = ndgrid(xx,yy); 

