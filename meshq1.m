function[E,x,y,nv,D,J] = meshq1(r,XX,YY,EE)
% for straight edged quadilateral ABCD
% X = 0.25*[ Va*(1-r)(1-s) + Vb*(1+r)(1-s) + Vc*(1+r)(1+s) + Va*(1-r)(1+s) ]
[R,S] = ndgrid(r); N = length(r); ne = length(EE(:,1));
R = reshape(R,[N*N,1]); S = reshape(S,[N*N,1]);
D = zeros(ne,N*N,4); J = zeros(ne,N*N); %D=[dxdr,dxds,dydr,dyds]
x = []; E = zeros(ne,N*N);
xx = zeros(N*N,2);
for ei = 1:ne
  KK = EE(ei,:);
  Vx = XX(KK);
  Vy = YY(KK);
  xx(:,1)=0.25*(Vx(1)*(1-R).*(1-S)+Vx(2)*(1+R).*(1-S)+Vx(3)*(1+R).*(1+S)+Vx(4)*(1-R).*(1+S));
  xx(:,2)=0.25*(Vy(1)*(1-R).*(1-S)+Vy(2)*(1+R).*(1-S)+Vy(3)*(1+R).*(1+S)+Vy(4)*(1-R).*(1+S));

  D(ei,:,1) = 0.25*(-Vx(1)*(1-S)+Vx(2)*(1-S)+Vx(3)*(1+S)-Vx(4)*(1+S));
  D(ei,:,2) = 0.25*(-Vx(1)*(1-R)-Vx(2)*(1-S)+Vx(3)*(1+R)+Vx(4)*(1-R));
  D(ei,:,3) = 0.25*(-Vy(1)*(1-S)+Vy(2)*(1-S)+Vy(3)*(1+S)-Vy(4)*(1+S));
  D(ei,:,4) = 0.25*(-Vy(1)*(1-R)-Vy(2)*(1-S)+Vy(3)*(1+R)+Vy(4)*(1-R));
  J(ei,:) = D(ei,:,1).*D(ei,:,4) - D(ei,:,2).*D(ei,:,3); J(ei,:) = abs(J(ei,:));
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