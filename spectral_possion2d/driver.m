%
lx1 = 10;
lxd = ceil(1.5*lx1);

[zm1,wm1] = zwgll(lx1-1);
[zmd,wmd] = zwgll(lxd-1); % dealias
[zme,wme] = zwgll(1);     % linear interpolation

J1d=interp_mat(zmd,zm1);
Je1=interp_mat(zm1,zme);
Jed=interp_mat(zmd,zme);

Dm1 = dhat(zm1);
Im1 = eye(lx1);
Dmd = dhat(zmd);
Imd = eye(lxd);

[xm1,ym1] = ndgrid(zm1);
[xmd,ymd] = ndgrid(zmd);

% deform geometry
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zm1,ze);
[xm1,ym1] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Je1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zmd,ze);
[xmd,ymd] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Jed);

% Jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Im1,Dm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Imd,Dmd);

% diag mass matrix
Bm1 = Jm1.*(wm1*wm1');
Bmd = Jmd.*(wmd*wmd');

% mask off dirichlet BC
msk = zeros(lx1,lx1);
msk(2:end-1,2:end-1) = 1;      % X: dir-dir, Y: dir-dir
msk(2:end-0,1:end-1) = 1;      % X: dir-neu, Y: neu-dir

% solve/ time step
f  = sin(pi*xm1).*sin(pi*ym1);
b  = mass(Bmd,Jd,msk,f);
ue = 0.5*f/pi/pi;

[Nx,Ny] = size(msk);
b = reshape(b,[Nx*Ny,1]);

% make conjugate gradient function

u = pcg(@(u)laplace(Dm1,Jd,Bmd,rxmd,rymd,sxmd,symd,msk,u),b,1e-8,100,[],[],b);

u = reshape(u,[Nx,Ny]);
% plt
mesh(xm1,ym1,u);

