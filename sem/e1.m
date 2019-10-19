%
lx1 = 10;
lxd = ceil(1.5*lx1);

[zm1,wm1] = zwgll(lx1-1);
[zmd,wmd] = zwgll(lxd-1); % dealias

Jd = interp_mat(zmd,zm1);

Dm1 = dhat(zm1);
Im1 = eye(lx1);
Dmd = dhat(zmd);
Imd = eye(lxd);

[xm1,ym1] = ndgrid(zm1);
[xmd,ymd] = ndgrid(zmd);

% deform geometry - gordon hall
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = g(lx1);
[xm1,ym1] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,lx1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = g(lxd);
[xmd,ymd] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,lxd);

% Jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac(xm1,ym1,Im1,Dm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac(xmd,ymd,Imd,Dmd);

% diag mass matrix
Bm1 = Jm1.*(wm1*wm1');
Bmd = Jmd.*(wmd*wmd');

% mask off dirichlet BC
msk = zeros(lx1,lx1);
msk(2:end-1,2:end-1) = 1;

% solve/ time step
f  = sin(pi*xm1).*sin(pi*ym1);
b = mass(Bmd,Jd,f);
ue = 0.5*f/pi/pi;

u = pcg(@(u)laplace(Dm1,Jd,Bmd,rxmd,rymd,sxmd,symd,u),b);

% plt


