%
% v <- kron(A,B) * u
% for A,B (m,n) matrices,
%       input  size n**2
%       output size m**2
%
function v = ABu(A,B,u)

	v = B*u*A';

end
