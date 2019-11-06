%
function [x,y] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,Ix,Iy,Jer,Jes)
%
% Jer - interpolate from e to r
% Jes - interpolate from e to s

xv = [xrm([1,end]),xrp([1,end])]; % vertices
yv = [yrm([1,end]),yrp([1,end])];

xv = ABu(Jes,Jer,xv);
yv = ABu(Jes,Jer,yv);

x = ABu(Jes,Ix,[xrm,xrp]) + ABu(Iy,Jer,[xsm';xsp']) - xv;
y = ABu(Jes,Ix,[yrm,yrp]) + ABu(Iy,Jer,[ysm';ysp']) - yv;

%=============================================================
if(0)
%------------------------------
figure;
fig=gcf;ax=gca;
hold on;grid on;
% title
title(['Mesh'],'fontsize',14);
% ax
ax.FontSize=14;
xlabel('$$X$$');
ylabel('$$Y$$');

mesh(x,y,0*x)
% color
colormap([0,0,0])
%------------------------------
figname=['mesh'];
saveas(fig,figname,'jpeg');
%------------------------------
end
%=============================================================

