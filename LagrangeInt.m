function[L] = LagrangeInt(xx)

% output Lagrange interpolating polynomails for interpolating points xx.
% N polynomials of order N-1.
N = length(xx);
D = diag(ones(N,1));
L = zeros(N);

x = linspace(xx(1),xx(end),1000); % fine mesh for plotting.
figure(); hold on; grid on;
for i = 1:N
  L(i,:) = polyfit(xx,D(:,i),N-1);
  y = polyval(L(i,:),x);
  plot(x,y,'linewidth',2);
end

end
