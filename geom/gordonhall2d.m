%
function [x,y] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zr,zs)
%
ze = [-1;1];
Jer = interp_mat(zr,ze);
Jes = interp_mat(zs,ze);
Ir  = eye(length(zr));
Is  = eye(length(zs));

xv = [xrm([1,end])';xrp([1;end])']; % vertices
yv = [yrm([1,end])';yrp([1;end])'];
xv = ABu(Jes,Jer,xv);
yv = ABu(Jes,Jer,yv);

x = ABu(Is,Jer,[xrm';xrp']) + ABu(Jes,Ir,[xsm,xsp]) - xv;
y = ABu(Is,Jer,[yrm';yrp']) + ABu(Jes,Ir,[ysm,ysp]) - yv;

%=============================================================
if(0)
%------------------------------
%figure;
fig=gcf;ax=gca;
hold off;grid on;
% title
title(['Mesh'],'fontsize',14);
% ax
ax.FontSize=14;
xlabel('$$X$$');
ylabel('$$Y$$');

mesh(x,y,0*x)
% color
colormap([0,0,0]);
view(2);
%------------------------------
figname=['mesh'];
saveas(fig,figname,'jpeg');
%------------------------------
end
%=============================================================

