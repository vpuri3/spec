% solve -(pu')' = f in 2 dimensions with arbitrary elements
clear all;
po=10;N=po+1;
addpath('../semhat'); [Bh,Dh,r,w] = semhat(po);

% make spectral mesh
%meshIllinois;
meshRect;
[E,x,y,Jac,GG,nv] = meshq1(r,XX,YY,EE); %Derev=[dxdr,dxds,dydr,dyds]

% input f, p, exact solution
f  = sin(pi*x).*sin(pi*y);
p  = ones(size(x));
ex = 0.5*f/(pi*pi);

% assemble
A = sparse(nv,nv);
B = sparse(nv,nv);
D = zeros(N*N,4);
G = zeros(N*N,3);
for ei = 1:ne
  K = E(ei,:);
  xe = x(K); ye = y(K); pe = p(K);
  J = Jac(ei,:);
  G(:,:) = GG(ei,:,:);
  %P = diag(pe); W = P*kron(Bh,Bh); I = eye(N);
  %Ae = kron(I,Dh') * W * kron(I,Dh) + kron(Dh',I) * W * kron(Dh,I);
  % rectangular
  P = diag(pe); W = P*kron(0.5*lx*Bh,0.5*ly*Bh); I = eye(N);
  Ae = kron(I,(2/lx)*Dh') * W * kron(I,(2/lx)*Dh) + kron((2/ly)*Dh',I) * W * kron((2/ly)*Dh,I);
  Be = (0.25*lx*ly)*kron(Bh,Bh);
  % Illinois
  %Ae = kron(I,(2/lx(ei))*Dh') * W * kron(I,(2/lx(ei))*Dh) + kron((2/ly(ei))*Dh',I) * W * kron((2/ly(ei))*Dh,I);
  %Be = kron(Bh,Bh);
  B(K,K) = B(K,K) + Be;
  A(K,K) = A(K,K) + Ae;
end

% hom dirichlet boundary restriction
% Illinois mesh
%tol=1e-12;
%Iflg = bitor(x==-1,x==2); Iflg = Iflg + bitor(y==0,y==5);
%Iflg = Iflg + bitand(y==1,x>=1-tol); Iflg = Iflg + bitand(y==1,x<=tol);
%Iflg = Iflg + bitand(y==4,x>=1-tol); Iflg = Iflg + bitand(y==4,x<=tol);
%Iflg = Iflg + bitand(y>=1-tol,x==1); Iflg = Iflg + bitand(y<=4+tol,x==1);
%Iflg = Iflg + bitand(y>=1-tol,x==0); Iflg = Iflg + bitand(y<=4+tol,x==0);
%Iflg = find(~Iflg); % indices corresponding to interior points.
% rectangular mesh
tol=1e-12;
Bflg = bitor(x==0,x==1);
Bflg = Bflg + bitor(y==0,y==1);
bflg = find(Bflg);  % boundary points
iflg = find(~Bflg); % interior points

% apply Dirichlet BC
I = speye(nv);
xb = x(bflg); yb = y(bflg);
szb = size(bflg); % # of dofs to remove
szi = size(iflg); % # of dofs to keep
ub = 0*xb;
A(bflg,iflg) = zeros(length(bflg),length(iflg)); A(bflg,bflg) = eye(length(bflg));
%B(bflg,iflg) = zeros(length(bflg),length(iflg));
h = B*f;
h(bflg) = ub;

% iterative solver setup
M = diag(diag(A));
maxit = 1e3; tol = 1e-13;
u = pcg(A,h,tol,maxit,M);
%u = A\h;

figure(); title('Numerical Solution u(x,y)'); hold on; addpath('../plt')
plot3(x,y,u,'ko');
grid on; shading interp; axis tight; colormap(fireice);

nx = sqrt(length(x));
X=reshape(x,[nx,nx]);
Y=reshape(y,[nx,nx]);
U=reshape(u,[nx,nx]);
%Er=reshape(u-ex,[nx,nx]);

figure(); title('Numerical Solution u(x,y)'); hold on; addpath('../plt')
surfc(X,Y,U);
grid on; shading interp; axis tight; colorbar; colormap(fireice);