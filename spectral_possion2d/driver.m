%
clear;
format compact; format shorte;

lx1 = 32;
lxd = ceil(1.5*lx1); lxd=lx1;

[zm1,wm1] = zwgll(lx1-1);
[zmd,wmd] = zwgll(lxd-1);    % dealias
[zme,wme] = zwgll(1);        % linear interpolation

J1d=interp_mat(zmd,zm1);
Je1=interp_mat(zm1,zme);
Jed=interp_mat(zmd,zme);

Dm1 = dhat(zm1);
Dmd = dhat(zmd);

Im1 = eye(lx1);
Imd = eye(lxd);

% deform geometry
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zm1,zme);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Je1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zmd,zme);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Jed);

%[xm1,ym1] = ndgrid(2*zm1);
%[xmd,ymd] = ndgrid(2*zmd);

% Jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Im1,Dm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Imd,Dmd);

% diag mass matrix
Bm1 = Jm1.*(wm1*wm1');
Bmd = Jmd.*(wmd*wmd');

% laplace operator
G11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
G12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
G22 = Bmd .* (sxmd.*sxmd + symd.*symd);

% mask off dirichlet BC
msk = zeros(lx1,lx1);
msk(2:end-1,2:end-1) = 1;      % X: dir-dir, Y: dir-dir
%msk(2:end-0,1:end-1) = 1;      % X: dir-neu, Y: neu-dir

% solve/ time step
f  = sin(pi*xm1).*sin(pi*ym1);
f  = 1+0*xm1;
b  = msk .* mass2d(Bmd,J1d,f);
ue = 0.5*f/pi/pi;

[u,iter,r2] = cg_possion2d(b,0*b,1e-8,500,Dm1,J1d,G11,G12,G22,msk);
mesh(xm1,ym1,u);
['iter:',num2str(iter),', res:',num2str(sqrt(r2)),' err:',num2str(max(max(abs(u-ue))))]

%[u] = laplace2d(f,Dm1,J1d,G11,G12,G22,msk);
%mesh(xm1,ym1,f); pause; mesh(xm1,ym1,u);

