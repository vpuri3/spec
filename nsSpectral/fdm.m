%
function [u] = fdm(b,Sr,Ss,Li);

	u = ABu(Ss',Sr',b);
	u = u .* Li;
	u = ABu(Ss,Sr,u);

end
