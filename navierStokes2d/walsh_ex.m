function [u,v,p] = walsh_ex(x,y,nu,t)
lam = 25;
u = exp(-  lam*nu*t)   * (sin(5*y)+     cos(3*x).*cos(4*y));
v = exp(-  lam*nu*t)   * (cos(5*x)+0.75*sin(3*x).*sin(4*y));   
p = exp(-2*lam*nu*t)/4 * (cos(7*x)+0.75*cos(7*y));

m = 4; n = 3; lam = m*m+n*n;
psi   =    sin(m*x).*cos(n*y);
psi_x =  m*cos(m*x).*cos(n*y);
psi_y = -n*sin(m*x).*sin(n*y);
u = exp(-  nu*lam*t)*(-psi_y);
v = exp(-  nu*lam*t)*( psi_x);
p = exp(-2*nu*lam*t)*0.5*(lam*psi.^2 + psi_x.^2 + psi_y.^2);

end
