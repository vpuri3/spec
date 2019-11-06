%
clear;
format compact; format shorte;
%----------------------------------------------------------------------
% todo
%
% inhomogeneous BC
% advection
%----------------------------------------------------------------------

lx1 = 8;
ly1 = 32;
lxd = ceil(1.5*lx1);
lyd = ceil(1.5*ly1);

[zrm1,wrm1] = zwgll(lx1-1);
[zsm1,wsm1] = zwgll(ly1-1);
[zrmd,wrmd] = zwgll(lxd-1);    % dealias
[zsmd,wsmd] = zwgll(lyd-1);
[zme,wme]   = zwgll(1);        % linear interpolation

Drm1 = dhat(zrm1);
Dsm1 = dhat(zsm1);
Drmd = dhat(zrmd);
Dsmd = dhat(zsmd);

Irm1 = eye(lx1);
Ism1 = eye(ly1);
Irmd = eye(lxd);
Ismd = eye(lyd);

Jr1d = interp_mat(zrmd,zrm1); % from lx1 to lxd
Js1d = interp_mat(zsmd,zsm1);
Jre1 = interp_mat(zrm1,zme);  % from e   to lx1
Jse1 = interp_mat(zsm1,zme);
Jred = interp_mat(zrmd,zme);  % from e   to lxd
Jsed = interp_mat(zsmd,zme);

%----------------------------------------------------------------------
% deform geometry
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zrm1,zsm1,zme);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Irm1,Ism1,Jre1,Jse1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zrmd,zsmd,zme);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Irmd,Ismd,Jred,Jsed);

% Jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);
%----------------------------------------------------------------------

% diag mass matrix
Bm1 = Jm1.*(wrm1*wsm1');
Bmd = Jmd.*(wrmd*wsmd');

% laplace operator
G11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
G12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
G22 = Bmd .* (sxmd.*sxmd + symd.*symd);

% mask off dirichlet BC
msk = zeros(lx1,ly1);
msk(2:end-1,2:end-1) = 1;      % X: dir-dir, Y: dir-dir
%msk(2:end-0,1:end-1) = 1;      % X: dir-neu, Y: neu-dir

%----------------------------------------------------------------------
% solution
f  = sin(pi*xm1).*sin(pi*ym1);
%f  = 1+0*xm1;
%ue = 0.5*f/pi/pi;

%----------------------------------------------------------------------
% RHS
b  = msk .* mass2d(Bmd,Jr1d,Js1d,f);

%----------------------------------------------------------------------
% solve

ifcg=0;
if(ifcg)
	[u,iter,r2] = cg_possion2d(b,0*b,1e-8,500,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22,msk);
else
	Rx = Irm1(2:end-1,:);
	Ry = Ism1(2:end-1,:);
	R  = sparse(kron(Ry ,Rx ));
	J  = sparse(kron(Js1d,Jr1d));
	
	G11=sparse(diag(reshape(G11,[lxd*lyd,1]))); G11 = J'*G11*J;
	G12=sparse(diag(reshape(G12,[lxd*lyd,1]))); G12 = J'*G12*J;
	G22=sparse(diag(reshape(G22,[lxd*lyd,1]))); G22 = J'*G22*J;
	G  = [G11 G12; G12 G22];
	
	Dr = kron(Ism1,Drm1);
	Ds = kron(Dsm1,Irm1);
	D  = [Dr;Ds];
	
	Ab = D'*G*D;
	A  = R*Ab*R';
	
	b  = R*reshape(b,[lx1*ly1,1]);
	
	u = A \ b;
	u = R'*u;
	u = reshape(u,[lx1,ly1]);
end
%----------------------------------------------------------------------
% plt
mesh(xm1,ym1,u);
