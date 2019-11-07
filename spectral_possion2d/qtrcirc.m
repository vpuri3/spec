%
function [xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = qtrcirc(zr,zs)
ze  = [-1;1];
Jer = interp_mat(zr,ze); % ze to zr
Jes = interp_mat(zs,ze);

xsm = zr*0.5;
ysm = zr*0 + 0.5;

asp = Jer*[3*pi/4;pi/4];
xsp  = cos(asp);
ysp  = sin(asp);

xrp = 0.5 + (zs+1)*(1/sqrt(2)-0.5)/2;
yrp = 0.5 + (zs+1)*(1/sqrt(2)-0.5)/2;

xrm = -xrp;
yrm =  yrp;

if(0)
	hold on; grid on;
	plot(xsm,ysm,'ro-','DisplayName','s-minus');
	plot(xsp,ysp,'rx-','DisplayName','s-plus ');
	plot(xrm,yrm,'bo-','DisplayName','r-minus');
	plot(xrp,yrp,'bx-','DisplayName','r-plus ');
	legend('show');
end
