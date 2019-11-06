%
function [xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = geom1(z,ze)
Je= interp_mat(z,ze);

xsm = z*0.5;
ysm = z*0 + 0.5;

asp = Je*[3*pi/4;pi/4];
xsp = cos(asp);
ysp = sin(asp);

xrp = 0.5 + (z+1)*(1/sqrt(2)-0.5)/2;
yrp = 0.5 + (z+1)*(1/sqrt(2)-0.5)/2;

xrm =-xrp;
yrm = yrp;

if(0)
	figure; hold on; grid on;
	plot(xsm,ysm,'ro-');      % s-minus
	plot(xsp,ysp,'mo-');      % s-plus
	plot(xrm,yrm,'bo-');      % r-minus
	plot(xrp,yrp,'ko-');      % r-plus
end
