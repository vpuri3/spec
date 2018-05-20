function[E,x,y,J,G,nv] = meshq1(r,XX,YY,EE)
% for parallelogram: V=0.25(Va+Vb+Vc+Vd)+(Vb-Va)*r+(Vd-Va)*s
% for any straight edged quadilateral ABCD
% EE(ei,:) is the list of vertices associated with element ei.
% XX(vi) contains X coordiate of vertex vi.

[R,S] = ndgrid(r); N = length(r); ne = length(EE(:,1));
R = reshape(R,[N*N,1]); S = reshape(S,[N*N,1]);
D = zeros(ne,N*N,4); % D=[dxdr,dxds,dydr,dyds]
J = zeros(ne,N*N);
G = zeros(ne,N*N,3); % G_11,G_12,G_22
x = []; E = zeros(ne,N*N);
xx = zeros(N*N,2);

for ei = 1:ne
  KK = EE(ei,:);
  Vx = XX(KK);
  Vy = YY(KK);
  xx(:,1)=0.25*sum(Vx)+0.5*(Vx(2)-Vx(1))*R+0.5*(Vx(4)-Vx(1))*S;
  xx(:,2)=0.25*sum(Vy)+0.5*(Vy(2)-Vy(1))*R+0.5*(Vy(4)-Vy(1))*S;
  % D=[dxdr,dxds,dydr,dyds]
  D(ei,:,1) = 0.5*(Vx(2)-Vx(1));
  D(ei,:,2) = 0.5*(Vx(4)-Vx(1));
  D(ei,:,3) = 0.5*(Vy(2)-Vy(1));
  D(ei,:,4) = 0.5*(Vy(4)-Vy(1));
  J(ei,:) = D(ei,:,1).*D(ei,:,4) - D(ei,:,2).*D(ei,:,3); J(ei,:) = abs(J(ei,:));
  G(ei,:,1) = D(ei,:,1).*D(ei,:,1) + D(ei,:,3).*D(ei,:,3);
  G(ei,:,2) = D(ei,:,1).*D(ei,:,2) + D(ei,:,3).*D(ei,:,4);
  G(ei,:,3) = D(ei,:,2).*D(ei,:,2) + D(ei,:,4).*D(ei,:,4);
  % follow lexicographical ordering within each element.
  if(ei==1) x = xx; end;
  x = [x; setdiff(xx,x,'rows')];
  %K = ismember(x,xx,'rows'); K = find(K);
  [~,K] = ismember(xx,x,'rows','R2012a');
  E(ei,:) = K;
end
nv = length(x);
y = x(:,2); x = x(:,1);


figure(); title('Quad Mesh'); hold on; grid on;
for ei = 1:ne
  K = E(ei,:);
  plot(x(K),y(K),'o','linewidth',0.8,'displayname',['Elem ' num2str(ei)]);
end
legend('show');

end