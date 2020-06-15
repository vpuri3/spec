%
% v <- kron(As,Br) * u
%
function v = ABu(As,Br,u)

	if(length(As)==0 & length(Br)==0); v=u;     return;
	elseif(length(As)==0);             v=Br*u ; return;
	elseif(length(Br)==0);             v=u*As'; return;
	end

	v = Br*u*As';

end
