% solve -(pu')' = f in 2 dimensions with 4 elements
po=20;N=po+1;
addpath('semhat'); [Bh,Dh,r,w] = semhat(po);
% rectangle mesh
%ne = 4; EE = zeros(ne,4); % vertices of the quad.
%EE(1,:) = [1,2,5,4]; EE(2,:) = [2,3,6,5]; EE(3,:) = [4,5,8,7]; EE(4,:) = [5,6,9,8];
%YY = [0,0,0,0.5,0.5,0.5,1,1,1]'; XX = [0,0.5,1,0,0.5,1,0,0.5,1];
% Illinois mesh
ne=7; EE=zeros(ne,4);
EE(1,:)=[1 2 15 16];
EE(2,:)=[2 3 6 15];
EE(3,:)=[3 4 5 6];
EE(4,:)=[15 6 7 14];
EE(5,:)=[13 14 11 12];
EE(6,:)=[14 7 10 11];
EE(7,:)=[7 8 9 10];
XX=[-1 0 1 2 2 1 1 2 2 1 0 -1 -1 0 0 -1];
YY=[0 0 0 0 1 1 4 4 5 5 5 5 4 4 1 1];
lx=ones(7,1);
ly=[1 1 1 2 1 1 1];
% make spectral mesh
[E,x,y,Derev,J] = meshq1(r,XX,YY,EE); %Derev=[dxdr,dxds,dydr,dyds]
nv = length(x);
% input f and p
f  = sin(pi*x).*sin(pi*y);
p  = ones(size(x));
ex = 0.5*f/pi/pi;

A = sparse(nv,nv);
B = sparse(nv,nv);
for ei = 1:ne
  K = E(ei,:);
  xe = x(K); ye = y(K); pe = p(K);
  P = diag(pe); W = P*kron(Bh,Bh); I = eye(N);
  %Ae = kron(I,Dh') * W * kron(I,Dh) + kron(Dh',I) * W * kron(Dh,I);
  Ae = kron(I,(2/lx(ei))*Dh') * W * kron(I,(2/lx(ei))*Dh) + kron((2/ly(ei))*Dh',I) * W * kron((2/ly(ei))*Dh,I);
  Be = kron(Bh,Bh);
  B(K,K) = B(K,K) + Be;
  A(K,K) = A(K,K) + Ae;
end
% hom dirichlet boundary restriction
tol=1e-12;
Iflg = bitor(x==-1,x==2); Iflg = Iflg + bitor(y==0,y==5);
Iflg = Iflg + bitand(y==1,x>=1-tol); Iflg = Iflg + bitand(y==1,x<=tol);
Iflg = Iflg + bitand(y==4,x>=1-tol); Iflg = Iflg + bitand(y==4,x<=tol);
Iflg = Iflg + bitand(y>=1-tol,x==1); Iflg = Iflg + bitand(y<=4+tol,x==1);
Iflg = Iflg + bitand(y>=1-tol,x==0); Iflg = Iflg + bitand(y<=4+tol,x==0);
Iflg = find(~Iflg); % indices corresponding to interior points.

% apply restrictions
A = A(Iflg,Iflg); B = B(Iflg,Iflg);
u = A \ B*f(Iflg);

ub = zeros(size(x)); ub(Iflg) = u; u = ub;
er = ex-u; max(abs(er))

figure(); title('Error u(x,y)'); hold on; addpath('~/Documents/MATLAB/plt')
plot3(x,y,u,'ko');
grid on; shading interp; axis tight; colormap(fireice);

X=reshape(x,[nx,nx]);
Y=reshape(y,[nx,nx]);
U=reshape(u,[nx,nx]);
Er=reshape(u-ex,[nx,nx]);

figure(); title('Error u(x,y)'); hold on; addpath('~/Documents/MATLAB/plt')
surfc(X,Y,Er);
grid on; shading interp; axis tight; colorbar; colormap(fireice);

figure(); title('Error u(x,y)'); hold on; addpath('~/Documents/MATLAB/plt')
surfc(X,Y,U);
grid on; shading interp; axis tight; colorbar; colormap(fireice);
