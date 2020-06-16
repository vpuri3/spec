%
function [xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = para(zr,zs)
ze  = [-1;1];
Jer = interp_mat(zr,ze); % ze to zr
Jes = interp_mat(zs,ze);

%-------------------
xsm = Jer*[-0.5; 0.5];
ysm = Jer*[-0.5;-0.5];
xsp = Jer*[-1.0; 1.0];
ysp = Jer*[ 0.5; 0.5];
ysp = Jer*[ 0.5; 0.5]; ysp = 1.0-0.5*xsp.*xsp;

xrm = Jer*[-0.5;-1.0];
yrm = Jer*[-0.5; 0.5];
xrp = Jer*[ 0.5; 1.0];
yrp = Jer*[-0.5; 0.5];

if(0)
	hold on; grid on;
	plot(xsm,ysm,'ro-','DisplayName','s-minus');
	plot(xsp,ysp,'rx-','DisplayName','s-plus ');
	plot(xrm,yrm,'bo-','DisplayName','r-minus');
	plot(xrp,yrp,'bx-','DisplayName','r-plus ');
	legend('show');pause;clf
end
