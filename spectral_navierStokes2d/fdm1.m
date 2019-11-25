%
function [u] = fdm1(b,RBi,Sr,Ss,Sri,Ssi,Li);

	u = b .* RBi;
	u = ABu(Ssi,Sri,u);
	u = u .* Li;
	u = ABu(Ss ,Sr ,u);

end
