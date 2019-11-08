%
clf;
format compact; format shorte;
%----------------------------------------------------------------------
%
% steady state advection-diffusion equation
%
%	nu*\del^2 u + q*u + \vect{c}\dot\grad{u}  = f
%   + Dirichlet/Neumann BC
%
% todo
%     -fast diagonalization
%
%----------------------------------------------------------------------

nx1 = 16;
ny1 = 16;
nxd = ceil(1.5*nx1);
nyd = ceil(1.5*ny1);

[zrm1,wrm1] = zwgll(nx1-1);
[zsm1,wsm1] = zwgll(ny1-1);
[zrmd,wrmd] = zwgll(nxd-1);
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

[xm1,ym1] = ndgrid(zrm1,zsm1);
[xmd,ymd] = ndgrid(zrmd,zsmd);
%----------------------------------------------------------------------
% data

solve=1; % 0: CG (possion), 1: FDM (possion), 2: direct (advection-diffusion)

nu= 1e-0;
q = 0+0*xm1;
cx= 0+0*xm1;
cy= 0+0*xm1;
f = sin(pi*xm1).*sin(pi*ym1);
f = 1+0*xm1;
%ue = 

% dirichlet BC
ub = 0-0*xm1;

% mask
Rx = Irm1(2:end-1,:);   % dir-dir
Ry = Ism1(2:end-1,:);   % neu-neu
msk = diag(Rx'*Rx) * diag(Ry'*Ry)';

ub = (1-msk) .* ub; % step function on dirichlet BC

%----------------------------------------------------------------------
% setup

% jacobian
[Jm1,Jim1,rxm1,rym1,sxm1,sym1] = jac2d(xm1,ym1,Irm1,Ism1,Drm1,Dsm1);
[Jmd,Jimd,rxmd,rymd,sxmd,symd] = jac2d(xmd,ymd,Irmd,Ismd,Drmd,Dsmd);

% mass
Bm1 = Jm1.*(wrm1*wsm1');
Bmd = Jmd.*(wrmd*wsmd');

% laplace operator setup
G11 = Bmd .* (rxmd.*rxmd + rymd.*rymd);
G12 = Bmd .* (rxmd.*sxmd + rymd.*symd);
G22 = Bmd .* (sxmd.*sxmd + symd.*symd);

% FDM - fast diagonalization setup
Lx = max(max(xm1))-min(min(xm1));
Ly = max(max(ym1))-min(min(ym1));

Br = (Lx/2)*diag(wrm1);
Bs = (Ly/2)*diag(wsm1);
Dr = (2/Lx)*Drm1;
Ds = (2/Ly)*Dsm1;
Ar = Dr'*Br*Dr;
As = Ds'*Bs*Ds;

Br = ABu(Rx,Rx,Br);
Bs = ABu(Ry,Ry,Bs);
Ar = ABu(Rx,Rx,Ar);
As = ABu(Ry,Ry,As);

[Sr,Lr] = eig(Ar,Br); Sri = inv(Sr);
[Ss,Ls] = eig(As,Bs); Ssi = inv(Ss);
%for j=1:length(Sr); Sr(:,j) = Sr(:,j)/sqrt(Sr(:,j)'*Sr(:,j)); end; % scaling err
%for j=1:length(Sr); Ss(:,j) = Ss(:,j)/sqrt(Ss(:,j)'*Ss(:,j)); end;
Lfdm = diag(Lr) + diag(Ls)';

%----------------------------------------------------------------------
% RHS
b = mass2d(Bmd,Jr1d,Js1d,f);							  % forcing
b = b - nu*laplace2d(ub,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22); % dirichlet bc
b = b - ub.*q;
b = b - advect2d(ub,cx,cy,Bmd,Irm1,Ism1,Jr1d,Js1d,Drm1,Dsm1,rxm1,rym1,sxm1,sym1);

b = msk .* b;
%----------------------------------------------------------------------
% solve

if(solve==0) % CG possion

	[uh,iter,r2] = cg_possion2d(b,0*b,1e-8,500,Jr1d,Js1d,Drm1,Dsm1,G11,G12,G22,msk);
	['cg iter:',num2str(iter),', residual:',num2str(sqrt(r2))]

elseif(solve==1) % FDM possion

	uh = (1/nu) * fdm(b,Bm1,Sr,Ss,Sri,Ssi,Rx,Ry,Lfdm);

elseif(solve==2) % direct advection-diffusion

	% data
	q  = reshape(q, [nx1*ny1,1]); Q = sparse(diag(q ));
	cx = reshape(cx,[nx1*ny1,1]);
	cy = reshape(cy,[nx1*ny1,1]);
	f  = reshape(f, [nx1*ny1,1]);
	b  = reshape(b, [nx1*ny1,1]); % rhs

	% operators
	R  = sparse(kron(Ry,Rx));
	J  = sparse(kron(Js1d,Jr1d));
	B  = sparse(diag(reshape(Bm1,[nx1*ny1,1])));
	Bd = sparse(diag(reshape(Bmd,[nxd*nyd,1])));

	% laplace op
	G11=sparse(diag(reshape(G11,[nxd*nyd,1]))); G11 = J'*G11*J;
	G12=sparse(diag(reshape(G12,[nxd*nyd,1]))); G12 = J'*G12*J;
	G22=sparse(diag(reshape(G22,[nxd*nyd,1]))); G22 = J'*G22*J;
	G  = [G11 G12; G12 G22];
	Dr = kron(Ism1,Drm1);
	Ds = kron(Dsm1,Irm1);
	D  = [Dr;Ds];
	A  = D'*G*D;

	% advection op
	RRX=sparse(diag(reshape(rxm1,[nx1*ny1,1])));
	RRY=sparse(diag(reshape(rym1,[nx1*ny1,1])));
	SSX=sparse(diag(reshape(sxm1,[nx1*ny1,1])));
	SSY=sparse(diag(reshape(sym1,[nx1*ny1,1])));
	Cxd=sparse(diag(J*cx));
	Cyd=sparse(diag(J*cy));
	Dx = RRX*Dr + SSX*Ds;
	Dy = RRY*Dr + SSY*Ds;
	C  = J'*Bd*(Cxd*J*Dx + Cyd*J*Dy);

	% system
	S  = R*(nu*A  + Q*B + C)*R';
	b  = R*b;
	uh = S \ b;
	uh = R'*uh;
	uh = reshape(uh,[nx1,ny1]);

end
%----------------------------------------------------------------------
% BC

u = uh + ub;

%----------------------------------------------------------------------
% plt
mesh(xm1,ym1,u);
xlabel('$$X$$');
ylabel('$$Y$$');
