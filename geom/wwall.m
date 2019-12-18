%
function [xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = wwall(zr,zs)
ze  = [-1;1];
Jer = interp_mat(zr,ze); % ze to zr
Jes = interp_mat(zs,ze);

xsm = Jer*[-0.5;0.5];
ysm = Jer*[-1  ;-1 ];

asp = Jer*[0;4*pi];
xsp = Jer*[-1;1];
ysp = 1+0.2*sin(asp);

xrp = Jes*[0.5;1];
yrp = Jes*[ -1;1];

xrm = -xrp;
yrm =  yrp;

if(0)
	hold on; grid on;
	plot(xsm,ysm,'ro-','DisplayName','s-minus');
	plot(xsp,ysp,'rx-','DisplayName','s-plus ');
	plot(xrm,yrm,'bo-','DisplayName','r-minus');
	plot(xrp,yrp,'bx-','DisplayName','r-plus ');
	legend('show');pause
end
