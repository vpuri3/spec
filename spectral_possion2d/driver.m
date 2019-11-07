%
clear;
format compact; format shorte;
%----------------------------------------------------------------------
% todo
%
% inhom, periodic bc
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

Drm1 = dhat(zrm1);
Dsm1 = dhat(zsm1);
Drmd = dhat(zrmd);
Dsmd = dhat(zsmd);

Irm1 = eye(lx1);
Ism1 = eye(ly1);
Irmd = eye(lxd);
Ismd = eye(lyd);

Jr1d = interp_mat(zrmd,zrm1); % lx1 to lxd
Js1d = interp_mat(zsmd,zsm1);

%----------------------------------------------------------------------
% deform geometry
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zrm1,zsm1);pause;
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);pause;

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zrmd,zsmd);
[xmd,ymd] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrmd,zsmd);

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

%----------------------------------------------------------------------
% solution
f  = sin(pi*xm1).*sin(pi*ym1);
f  = 1+0*xm1;
%ue = 

% mask off dirichlet BC
msk = zeros(lx1,ly1);
msk(2:end-1,2:end-1) = 1;      % X: dir-dir, Y: dir-dir
msk(1:end-1,2:end-1) = 1;      % X: neu-dir, Y: dir-dir
%msk(2:end-0,1:end-1) = 1;      % X: dir-neu, Y: neu-dir

%----------------------------------------------------------------------
% RHS
b  = msk .* mass2d(Bmd,Jr1d,Js1d,f);

%----------------------------------------------------------------------
% solve

ifcg=1;
if(ifcg)
	[u,iter,r2] = cg_possion2d(b,0*b,1e-8,500,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22,msk);
	['cg iter:',num2str(iter),', residual:',num2str(sqrt(r2))]
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
