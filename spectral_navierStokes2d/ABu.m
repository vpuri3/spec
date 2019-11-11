%
% v <- kron(As,Br) * u
%
function v = ABu(As,Br,u)

	v = Br*u*As';

end
