% solve u_xx = exp(u) with hom. Dir. BC using fixed point iteration
clf, hold on; grid on;

N = 50; k = 9;
[D,xc] = cheb(N);
D2 = D*D; D2 = D2(2:N,2:N);
I = eye(size(D2));
[xx,yy] = meshgrid(xc(2:N));
x = xx(:); y = yy(:);
f = 10*sin(8*x.*(y-1));
L = kron(I,D2) + kron(D2,I) + k*k*kron(I,I);

figure(1), clf, spy(L);
tic, u = L\f; toc

u = reshape(u,N-1,N-1);
uu = zeros(N+1); uu(2:N,2:N) = u;
[xx,yy] = meshgrid(xc);
[xxx,yyy] = meshgrid(linspace(-1,1,100));
uuu = interp2(xx,yy,uu,xxx,yyy);

figure(2), clf, mesh(xxx,yyy,uuu), colormap([0 0 0]);
xlabel x, ylabel y, zlabel u;

figure(3), clf, contour(xxx,yyy,uuu), colormap([0 0 0]), axis square
xlabel x, ylabel y;
