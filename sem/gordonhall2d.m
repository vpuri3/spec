%
function [x,y] =  gordonhall2d(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,lx)
%
% Godron Hall
% refer to
% https://wiki.math.ntnu.no/_media/ma8502/2014h/deformed2.pdf
%

[z ,w ] = zwgll(lx-1);
[ze,we] = zwgll(1);

I = eye(lx);
Ie= eye(2);
Je= interp_mat(z,ze);

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

