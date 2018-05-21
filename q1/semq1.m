% solve -(pu')' = f in 2 dimensions with arbitrary elements
clear all;
po=20;N=po+1;
addpath('../semhat'); [Bh,Dh,r,w] = semhat(po);

% make GLL element mesh
%meshIllinois;
meshRect;
[E,x,y,J,G,nv] = meshq1(r,XX,YY,EE);

% input f, p, exact solution
f  = sin(pi*x).*sin(pi*y);
p  = ones(size(x));
ex = 0.5*f/(pi*pi);

% assemble
A = sparse(nv,nv);
B = sparse(nv,nv);
Je = zeros(N*N,1);
Ge = zeros(N*N,3);
% assemble and turn into an iterative solver. make function afun
W = kron(w,w); I = eye(N);
for ei = 1:ne
  K = E(ei,:);
  xe = x(K); ye = y(K); pe = p(K);
  Je(:,:) = J(ei,:);
  Ge(:,:) = G(ei,:,:);
  Ge(:,1) = pe.*Ge(:,1); Ge(:,2) = pe.*Ge(:,2); Ge(:,3) = pe.*Ge(:,3);
  Ae =    kron(I,Dh')*diag(W.*Ge(:,1))*kron(I,Dh);
  Ae = Ae+kron(I,Dh')*diag(W.*Ge(:,2))*kron(Dh,I);
  Ae = Ae+kron(Dh',I)*diag(W.*Ge(:,2))*kron(I,Dh);
  Ae = Ae+kron(Dh',I)*diag(W.*Ge(:,3))*kron(Dh,I);
  %rectangular
  %P = diag(pe); W = P*kron(0.5*lx*Bh,0.5*ly*Bh)
  %Ae = kron(I,(2/lx)*Dh') * W * kron(I,(2/lx)*Dh) + kron((2/ly)*Dh',I) * W * kron((2/ly)*Dh,I);

  Be = diag(Je)*kron(Bh,Bh);
  B(K,K) = B(K,K) + Be;
  A(K,K) = A(K,K) + Ae;
end

% hom dirichlet boundary condition
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
er = max(abs(ex-u));

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
