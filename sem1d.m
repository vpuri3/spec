% solve -(pu')' = f
addpath('semhat');
p=16;N=p+1; % # of points per element
[Bh,Dh,r,w] = semhat(p);
P = eye(N);
Ah = Dh'*(P*Bh)*Dh;

ne = 2;
L = 1; l = L/ne;
xx = linspace(0,L,ne+1); h = diff(xx);

x = 0.5*( xx(2)+xx(1)+(xx(2)-xx(1))*r ); h = [1];
for i = 3:length(xx)
  a = 0.5*( xx(i-1)+xx(i)+(xx(i)-xx(i-1))*r );
  x = [x(1:end-1);a];
  h = [h; zeros(N-2,1); 1];
end

%%%%%% element to global
A = (2/l)*conv2(diag(h),Ah); % global neumann operator.
B = (l/2)*conv2(diag(h),Bh);
f = exp(x);

n = length(x); % # of global degrees of freedom.
AA = zeros(ne*N); Q = zeros(ne*N,n);
AA(1:N,1:N) = Ah; Q(1:N,1:N) = eye(N);
for i = 1:(ne-1)
  AA(i*N+(1:N),i*N+(1:N)) = Ah;
  Q(i*N+(1:N),i*N-i+(1:N)) = eye(N);
end
AA = (2/l)*AA; % A == Q' * A * Q ; global Neumann operator.


%%%%% apply restrictions
% R = [0 I 0]; u_dirichlet = R * u_global;
% u_global = R' * u_dirichlet
% A_dirichlet = R * A_global * R'
% B_dirichlet = R * B
dirichlet = 1;
if(dirichlet) % hom dirichlet BC on right edge
  A = A(2:end-1,2:end-1); B = B(2:end-1,2:end-1);
  u = A\B*f(2:end-1); u=[0;u;0];
end
if(~dirichlet) % restrict solution such that average is zero.
  A(end,:) = A(end,:) + ones(1,n);
  u = A \ B*f;
end





figure(); title('Solution')
hold on; grid on;
plot(x,u,'-','linewidth',2,'displayname','u');
plot(x,zeros(size(x)),'o','linewidth',1.2,'displayname','Points');
legend('show')


% need some sort of a numbering to create a function.
% solve poisson's equation on a rectangle.