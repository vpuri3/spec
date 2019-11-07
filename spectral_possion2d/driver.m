%
clf;
format compact; format shorte;
%----------------------------------------------------------------------
%
% solves -\del^2 u = f
%        + Dirichlet/Neumann BC
%
%----------------------------------------------------------------------

nx1 = 8;
ny1 = 16;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);

[zrm1,wrm1] = zwgll(nx1-1);
[zsm1,wsm1] = zwgll(ny1-1);
[zrmd,wrmd] = zwgll(nxd-1);    % dealias
[zsmd,wsmd] = zwgll(nyd-1);

Drm1 = dhat(zrm1);
Dsm1 = dhat(zsm1);
Drmd = dhat(zrmd);
Dsmd = dhat(zsmd);

Irm1 = eye(nx1);
Ism1 = eye(ny1);
Irmd = eye(nxd);
Ismd = eye(nyd);

Jr1d = interp_mat(zrmd,zrm1); % nx1 to nxd
Js1d = interp_mat(zsmd,zsm1);

%----------------------------------------------------------------------
% deform geometry
[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrm1,zsm1);
[xm1,ym1] = gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zrm1,zsm1);

[xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zrmd,zsmd);
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
% data
f  = sin(pi*xm1).*sin(pi*ym1);
f  = 1+0*xm1;
%ue = 

% mask off dirichlet BC
msk = zeros(nx1,ny1);
msk(2:end-1,2:end-0) = 1;      % X: neu-dir, Y: dir-dir

% dirichlet BC
ub = 0+0.05*xm1;
ub = (1-msk) .* ub;

%----------------------------------------------------------------------
% RHS
b  = mass2d(Bmd,Jr1d,Js1d,f);
b  = b - laplace2d(ub,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22);

b  = msk .* b;
%----------------------------------------------------------------------
% solve

ifcg=0;
if(ifcg)
	[uh,iter,r2] = cg_possion2d(b,0*b,1e-8,500,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22,msk);
	['cg iter:',num2str(iter),', residual:',num2str(sqrt(r2))]
else
	Rx = Irm1(2:end-1,:);
	Ry = Ism1(2:end-0,:);
	R  = sparse(kron(Ry,Rx));
	J  = sparse(kron(Js1d,Jr1d));
	
	G11=sparse(diag(reshape(G11,[nxd*nyd,1]))); G11 = J'*G11*J;
	G12=sparse(diag(reshape(G12,[nxd*nyd,1]))); G12 = J'*G12*J;
	G22=sparse(diag(reshape(G22,[nxd*nyd,1]))); G22 = J'*G22*J;
	G  = [G11 G12; G12 G22];
	
	Dr = kron(Ism1,Drm1);
	Ds = kron(Dsm1,Irm1);
	D  = [Dr;Ds];
	
	Ab = D'*G*D;
	A  = R*Ab*R';
	
	b  = R*reshape(b,[nx1*ny1,1]);
	
	uh = A \ b;
	uh = R'*uh;
	uh = reshape(uh,[nx1,ny1]);
end
%----------------------------------------------------------------------
% BC

u = uh + ub;

%----------------------------------------------------------------------
% plt
mesh(xm1,ym1,u);
