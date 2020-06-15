%
function x = geomspace(x0,x1,n,s)

% linspace but with geom progression (scaling factor s)

if(s==1.0) x=linspace(0,1,n); return; end;

x=zeros(n,1);
x(2)=(s-1)/(s^(n-1)-1);

for i=3:n
	dx  =x(i-1)-x(i-2);
	x(i)=x(i-1)+dx*s;
end

x=x0+x*(x1-x0);

%clf;plot(x,0*x,'ko'); grid on;
