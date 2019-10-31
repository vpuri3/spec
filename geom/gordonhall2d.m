%
function [x,y] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Je)
%
% Je - interpolation matrix from {-1,1} to basis points of length lx
%

lx = size(Je,1);
I = eye(lx);

xv = [xrm([1,end]),xrp([1,end])]; % vertices
yv = [yrm([1,end]),yrp([1,end])];

xv = ABu(Je,Je,xv);
yv = ABu(Je,Je,yv);

x = ABu(Je,I,[xrm,xrp]) + ABu(I,Je,[xsm';xsp']) - xv;
y = ABu(Je,I,[yrm,yrp]) + ABu(I,Je,[ysm';ysp']) - yv;

%=============================================================
if(0)
%------------------------------
figure;
fig=gcf;ax=gca;
hold on;grid on;
% title
title(['',' Mesh',' lx = ',num2str(lx)],'fontsize',14);
% ax
ax.FontSize=14;
xlabel('$$X$$');
ylabel('$$Y$$');

mesh(x,y,0*x)
% color
colormap([0,0,0])
%------------------------------
figname=['','mesh','lx',num2str(lx)];
saveas(fig,figname,'jpeg');
%------------------------------
end
%=============================================================

