function [a,b] = bdfext3(t)
%
%	t - array of time-steps (most recent to oldest)
%
%	a(3  ,1)
%	b(3+1,1)
%
%
	[a,b] = bdfextk(t);
    k = length(a);     % order

	if(k<3)

		a = [a;zeros(3-k,1)];
		b = [b;zeros(3-k,1)];

	else

		a = a(1:3);
		b = b(1:4);

	end

end

