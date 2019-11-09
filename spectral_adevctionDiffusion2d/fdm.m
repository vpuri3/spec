%
function [u] = fdm(b,Bi,Sr,Ss,Sri,Ssi,Rx,Ry,Di);

	u = b .* Bi;
	u = ABu(Ry ,Rx ,u);
	u = ABu(Ssi,Sri,u);
	u = u .* Di;
	u = ABu(Ss ,Sr ,u);
	u = ABu(Ry',Rx',u);

end
