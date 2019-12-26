function [ue,ve] = kov_ex(x,y,Re)

lam = Re/2 - sqrt(Re^2/4 + 4*pi^2);
ue = 1 - exp(lam.*x).*cos(2*pi.*y);
ve = lam/(2*pi) * exp(lam*x).*sin(2*pi*y);

end
