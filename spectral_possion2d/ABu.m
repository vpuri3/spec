%
% v <- kron(As,Br) * u
% for As,Br (m,n) matrices,
%       input  size n**2
%       output size m**2
%
function v = ABu(As,Br,u)

	v = Br*u*As';

end
