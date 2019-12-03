% Geometry
a = 0; Lx = 2*pi; Ly=2*pi;
x = a + Lx/2 * (zm1+1); y = a + Ly/2 * (zm1+1); [xm1,ym1]=ndgrid(x,y);
x = a + Lx/2 * (zm2+1); y = a + Ly/2 * (zm2+1); [xm2,ym2]=ndgrid(x,y);
x = a + Lx/2 * (zmd+1); y = a + Ly/2 * (zmd+1); [xmd,ymd]=ndgrid(x,y);

%initial conditions
[Ux,Uy,p0] = walsh_ex(X,Y,nu,0); %t=0, check pressure exact sol in walsh_ex!!

% Dir-Dir BC
Rx=Rx(2:end-1,:);
Ry=Ry(2:end-1,:);

%time dependent BC

[vxe,vye,~] = walsh_ex(xm1,ym1,visc1,t);
%left face: x=0
vxb(1,:) = vxe(1,:);
vyb(1,:) = vye(1,:);
%right face: x=2*pi %can be periodic
vxb(end,:) = vxe(end,:);
vyb(end,:) = vye(end,:);
%bottom face: y=0
vxb(:,1) = vxe(:,1);
vyb(:,1) = vye(:,1);
%top face: y=2*pi %can be periodic
vxb(:,end) = vxe(:,end);
vyb(:,end) = vye(:,end);
