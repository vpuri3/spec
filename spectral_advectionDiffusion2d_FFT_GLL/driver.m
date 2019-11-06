function [uegy,uinf,time] = driver(Nx,Ny,T,ifadv);
%
% [ued,uid,td] = hw(60,20,1e2,0);
% [uea,uia,ta] = hw(60,20,1e2,1);
%
format compact; format shorte;

% fourier
lx=2*pi; hx=lx/Nx; x=hx*[0:Nx-1]';
Ix=speye(Nx); Rx=Ix; nx=Nx;
Dx=(2*pi/lx)*dhatf(Nx); Dhh=Dx*Dx;  
Bx=pi*Ix; Bx(1,1)=2*pi;  Bx=lx*Bx/(2*pi); wx=diag(Bx);
Ax=-Bx*Dhh; 

% legendre
ly=2;
[By,Dy,y,wy] = semhat(Ny);   Ay=Dy'*By*Dy;
Iy=speye(Ny+1);              Ry=Iy;Ry=Ry(2:end-1,:); ny=size(Ry,1);
By=(ly/2)*Ry*By*Ry';         Ay=(2/ly)*(Ry*Ay*Ry');

% dealias
Mx=ceil(1.5*Nx);
My=ceil(1.5*Ny);
hxd=lx/Mx;xd=hxd*[0:Mx-1]';
wxd=pi*ones(Mx,1);wxd(1)=2*pi; wxd=wxd*lx/(2*pi);
[yd,wyd]=zwgll(My);
Jx=interp_mat(xd,x);
Jy=interp_mat(yd,y);
[Xd,Yd]=ndgrid(xd,yd);

% diagonal mass matrices
B =wx*wy';
RB=Rx*B*Ry';
Bd=wxd*wyd';

% fast diagonalization setup
[Sx,Lx]=eig(full(Ax),full(Bx));ex=ones(nx,1);
[Sy,Ly]=eig(full(Ay),full(By));ey=ones(ny,1);
for j=1:ny; Sy(:,j)=Sy(:,j)/sqrt(Sy(:,j)'*By*Sy(:,j)); end;

% initialize
[X,Y]=ndgrid(x,y);
Cx=0*Yd; if(ifadv);Cx=1-Yd.*Yd;end;  % advection speed X
Cy=0*Xd;                             % advection speed Y
nu=0.001;                            % diffusivity

% time stepper
dt=1e-2;
nt=T/dt;

% initial condition
t=0;
ub=min(X>1,X<2);
u = Rx*ub*Ry';
uh=rfft(u);

uh0=uh*0;
uh1=uh0;
uh2=uh0;
uh3=uh0;
f1 =uh0;
f2 =uh0;

% norms
uinf = zeros(nt,1);
uegy = zeros(nt,1);
time = zeros(nt,1);

% test pure advection
%
%Cx=0*Yd+1;Cy=0*Xd+0; nu=0;
%ub=sin(X).*sin(pi*Y); uh=Rx*rfft(ub)*Ry'; ue=ub;
%ub=min(X>1,X<2);      uh=Rx*rfft(ub)*Ry'; ue=ub;

% test steady diffusion
%
%b=0*b;nu=1;
%scale = 1+pi*pi;Ue=sin(X).*sin(pi*Y); F=scale*Ue;
%r=Rx*F*Ry'; r=RB.*rfft(r);

% advection no aliasing
%f1 = Rx'*uh*Ry;   f1 = Cx.*irfft(Dx*f1) + Cy.*irfft(f1*Dy');
%f1 = B.*rfft(f1); f1 = -Rx*f1*Ry';

for i=1:nt
	t = t+dt;
	[a,b] = bdfex3(i,dt);
	
	uh3=uh2;uh2=uh1;uh1=uh;
	f3=f2;f2=f1;
	f1= -Rx*advec_op(Rx'*uh*Ry,Cx,Cy,Bd,Dx,Dy,Jx,Jy)*Ry';

	r = a(1)*f1 +a(2)*f2 +a(3)*f3 - RB.*(b(2)*uh1+b(3)*uh2+b(4)*uh3);

	if(i<=3);
		D= b(1) + nu*( diag(Lx)*ey' + ex*diag(Ly)' );
		D=1./D;
	end

	% solve
	uh = (1./RB).*r;
	uh = Sx'*r*Sy;
	uh = uh .* D;
	uh = Sx*uh*Sy';

	% norm
	ubh= Rx'*uh*Ry;
	ub = irfft(ubh);
	uinf(i) = max(max(abs(ub)));
	uegy(i) = hx*sum(sqrt(ub.*ub)*wy) / (lx*ly);
	%uegy(i) = wx'*(sqrt(ubh.*ubh)*wy) / (lx*ly);
	time(i) = t;

	% plotting
	if(mod(i,0.1*nt)==0);
		hold off;
		mesh(X,Y,ub);title(['t=',num2str(t)]);pause(0.05);
	end

end
%max(max(abs(ue-ub)))

%=============================================================
if(0)
%------------------------------
figure;
fig=gcf;ax=gca;
hold on;grid on;
title(['Diffusion $$T=$$',num2str(t)],'fontsize',14);
%lgd=legend('location','northwest');lgd.FontSize=10;
% ax
ax.XScale='linear';
ax.YScale='linear';
xlabel('$$x$$');
ylabel('$$y$$');

surf(X,Y,ub); shading interp; view(3);
%------------------------------
figname=['d','_','t_',num2str(t)];
saveas(fig,figname,'jpeg');
%------------------------------
end
%=============================================================
if(0)
%------------------------------
figure;
fig=gcf;ax=gca;
hold on;grid on;
% title
title(['Diffusion $$\|u\|_2,\|u\|_\infty$$'],'fontsize',14);
%lg2
lgd=legend('location','northeast');lgd.FontSize=10;
% ax
ax.XScale='linear';
ax.YScale='linear';
xlabel('$$t$$');
ylabel('$$\|u\|_2,\|u\|_\infty$$');

plot(time,uinf,'-','linewidth',2.0,'DisplayName','$$\|u\|_\infty$$');
plot(time,uegy,'-','linewidth',2.0,'DisplayName','$$\|u\|_2     $$');
%------------------------------
figname=['d','_','norm'];
saveas(fig,figname,'jpeg');
%------------------------------
end
%=============================================================
