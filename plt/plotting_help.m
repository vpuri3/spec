% Plotting
close all;
hold on % Plot on the same figure

%---------------------------------------------------------------
plot(x,y,'bo',var,P,'r--')  % Plot...add as many x and y combinations as needed 
legend('Data Points','Langrange Approximation') % names in order of plotting above
% Colors: b = blue, k = black, r = red, g = green, y = yellow
% Styling: o = circles, - = solid line, -- = dashed line


loglog , semilogx , semilogy % logarithmic, semi-logarithmic


% Plotting with loops eg:
hold on
for c = 1:d
    len = 5/h(c)+1;
    plot(T(1:len,c),er_F(1:len,c),'-o','Linewidth',1.75,'DisplayName',['h=' num2str(h(c))])
end
title('Absolute Error in Forward Euler with varied grid spacing (h)','FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Absolute Error','FontSize',14)
axis sqaure
legend('show')
%---------------------------------------------------------------
% plotting
figure; fig=gcf; ax=gca;
ax.XScale='log';
ax.YScale='linear';
ax.FontSize=14;
xlabel('y-plus')
ylabel('u-plus')

title(['SmoothWavyWall, atime: ' num2str(atime)],'fontsize',14)
hold on; grid on; axis square;

semilogx(yp(V),up(V),'ko','linewidth',2.0,'displayname','Valley');
semilogx(yp(P),up(P),'ro','linewidth',2.0,'displayname','Peaks');

legend('location','northwest')
%---------------------------------------------------------------

%=============================================================
if(1)
%------------------------------
figure;
fig=gcf;ax=gca;
hold on;grid on;
%clf(fig);colordef(fig,'black');
% title
title([casename,' TKE Turbulent Diffusion',' atime:',num2str(at),'s'],'fontsize',14);
% pos
daspect([1,1,1]);set(fig,'position',[0,0,1000,500])
% ax
xlabel('x/\lambda_1');
ylabel('y');
xlim([min(min(x)),max(max(x))]);
ylim([min(min(y)),max(max(y))]);

pcolor(x,y,tdK);
% color
shading interp;colormap jet;
hcb=colorbar;title(hcb,['\eta_{ij}/U^2'],'fontsize',14);
if(cb) caxis([budgetsmin,budgetsmax]);end;
%------------------------------
figname=[cname,'-','tke-td'];
saveas(fig,figname,'jpeg');
end
%=============================================================
if(1)
%------------------------------
figure;
fig=gcf;ax=gca;
hold on;grid on;
% title
title(['Advection Diffusion Convergence'],'fontsize',14);
%lg2
lgd=legend('location','southwest');lgd.FontSize=10;
% ax
ax.XScale='log';
ax.YScale='log';
xlabel('n');
ylabel('error');

for j=1:length(a)
	p=plot(n,er(:,j),'-','linewidth',2.0);
	p.DisplayName=['Pe = ',num2str(Pe(j))];
end
p=plot(order,order_2,'k--+','linewidth',2.0);
p.DisplayName=['Second Order Line'];
%------------------------------
figname=['conv','-','ad'];
saveas(fig,figname,'jpeg');
%------------------------------
end
%=============================================================
