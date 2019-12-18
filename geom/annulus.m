%
function [xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp] = annulus(r0,r1,span,zr,zs)
ze  = [-1;1];
Jer = interp_mat(zr,ze); % ze to zr
Jes = interp_mat(zs,ze);

xrp = 0*zs;
yrp = Jes*[-r0;-r1];

spn = span-pi/2;

xrm = Jes*[r0;r1]*cos(spn);
yrm = Jes*[r0;r1]*sin(spn);

as = Jer*[-pi/2;spn];
xsm = r0*cos(as);
ysm = r0*sin(as);

xsp = r1*cos(as);
ysp = r1*sin(as);

if(0)
	hold on; grid on;
	plot(xsm,ysm,'ro-','DisplayName','s-minus');
	plot(xsp,ysp,'rx-','DisplayName','s-plus ');
	plot(xrm,yrm,'bo-','DisplayName','r-minus');
	plot(xrp,yrp,'bx-','DisplayName','r-plus ');
	legend('show');pause
end
