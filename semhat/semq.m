
function[Q] =  semq(E,N,bc)

m=N+1; Q=spalloc(E*m,E*N+1,3*E*m);

i=1; j=1;
for e=1:E; Q(i:(i+N),j:(j+N))=speye(m); i=i+m; j=j+N; end;

if bc==2; Q(end,1)=1; Q=Q(:,1:end-1); end;  %% Periodic BC

