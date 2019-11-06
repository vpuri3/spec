%
function [xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(zr,zs,ze)
Jer = interp_mat(zr,ze); % from ze to zx
Jes = interp_mat(zs,ze);

xsm = zs*0.5;
ysm = zs*0 + 0.5;

asp = Jes*[3*pi/4;pi/4];
xsp  = cos(asp);
ysp  = sin(asp);

xrp = 0.5 + (zr+1)*(1/sqrt(2)-0.5)/2;
yrp = 0.5 + (zr+1)*(1/sqrt(2)-0.5)/2;

xrm = -xrp;
yrm =  yrp;

if(0)
	figure; hold on; grid on;
	plot(xsm,ysm,'ro-','DisplayName','s-minus');
	plot(xsp,ysp,'mo-','DisplayName','s-plus ');
	plot(xrm,yrm,'bo-','DisplayName','r-minus');
	plot(xrp,yrp,'ko-','DisplayName','r-plus ');
	legend('show');
end
