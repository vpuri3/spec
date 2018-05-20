% solve -(pu')' = f in 2D with rectangular elements.
po=16;N=po+1; nex=3; ney=3; ne=nex*ney;
addpath('semhat'); [Bh,Dh,r,w] = semhat(po);
L = 1; lx = L/nex; ly = L/ney;

xx = linspace(0,1,nex+1); yy = linspace(0,1,ney+1); % grid demarkations
[x,y,E,nv]=meshrect(r,xx,yy);

f  = sin(pi*x).*sin(pi*y);
p  = ones(size(x));
ex = 0.5*f/pi/pi;

% assemble into sparse matrix.
A = sparse(nv,nv);
B = sparse(nv,nv);
for ei = 1:ne
  K = E(ei,:);
  xe = x(K); ye = y(K); pe = p(K);
  P = diag(pe); W = P*kron(0.5*lx*Bh,0.5*ly*Bh); I = eye(N);
  Ae = kron(I,(2/lx)*Dh') * W * kron(I,(2/lx)*Dh) + kron((2/ly)*Dh',I) * W * kron((2/ly)*Dh,I);
  Be = (0.25*lx*ly)*kron(Bh,Bh);
  B(K,K) = B(K,K) + Be;
  A(K,K) = A(K,K) + Ae;
end

% hom dirichlet boundary restriction
Iflg = bitor(x==0,x==1); Iflg = Iflg + bitor(y==0,y==1);
Iflg = find(~Iflg); % indices corresponding to interior points.

% apply restrictions
A = A(Iflg,Iflg); B = B(Iflg,Iflg);
u = A \ B*f(Iflg);

ub = zeros(size(x)); ub(Iflg) = u; u = ub;
er = ex-u; max(abs(er))

nx = sqrt(nv);
X=reshape(x,[nx,nx]);
Y=reshape(y,[nx,nx]);
U=reshape(u,[nx,nx]);
Er=reshape(u-ex,[nx,nx]);

figure(); title('Solution u(x,y)'); hold on; addpath('~/Documents/MATLAB/plt')
surfc(X,Y,U);
grid on; shading interp; axis tight; colorbar; colormap(fireice);

figure(); title('Error'); hold on; addpath('~/Documents/MATLAB/plt')
surfc(X,Y,Er);
grid on; shading interp; axis tight; colorbar; colormap(fireice);