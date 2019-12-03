function [a,b] = bdfextk(t)
%
%	t - array of time-steps (most recent to oldest)
%
%	a(k  ,1)
%	b(k+1,1)
%
%
	t = unique(t,'stable');
    k = length(t)-1;

	tcurr = t(1);
	tprev = t(2:end);
	a = interp_mat(tcurr,tprev)'; % extrapolate

	D = deriv_mat(t);
	b = D(1,:)';
