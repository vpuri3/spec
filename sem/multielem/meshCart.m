% constructs rectangular Spectral Element mesh.
% INPUT: points in canonical domain [-1,1], grid demarcations in x and y
% direction
% ASSUMES z contains end points.
function [x,y,E,nv] = meshrect(z,xx,yy)
nex = length(xx)-1; ney = length(yy)-1;
ne = nex*ney; N = length(z);
gx = []; gy = [];
for i = 1:max(nex,ney)
  j = i+1;
  if(i<=nex)
    ax = 0.5*(xx(j-1)+xx(j)+(xx(j)-xx(j-1))*z);
    gx = [gx(1:end-1);ax];
  end
  if(i<=ney)
    ay = 0.5*(yy(j-1)+yy(j)+(yy(j)-yy(j-1))*z);
    gy = [gy(1:end-1);ay];
  end
end
[X,Y] = ndgrid(gx,gy); nv = length(gx)*length(gy);
x = reshape(X,[nv,1]); y = reshape(Y,[nv,1]);
E = zeros(ne,N*N); % lexicographical ordering of elements
tol = 1e-12;

figure(); title('Rectangle Mesh'); hold on; grid on;
for j = 1:ney
for i = 1:nex
  ei = i+(j-1)*nex;
  s = find(x>=xx(i)-tol & x<=xx(i+1)+tol);
  t = find(y>=yy(j)-tol & y<=yy(j+1)+tol);
  K = intersect(s,t)';
  E(ei,:) = K;
  plot(x(K),y(K),'o','linewidth',0.8,'displayname',['Elem ' num2str(ei)]);
end
end
legend('show');
end
