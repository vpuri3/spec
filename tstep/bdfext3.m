function [a,b] = bdfext3(t)
%
%	t - sequence of time-steps (lastest to oldest)
%
%	a(k  ,1)
%	b(k+1,1)
%
%
    k = length(t)-1;

	tcurr = t(1);
	tprev = t(2:end);
	a = interp_mat(tcurr,tprev)';

	D = deriv_mat(t);
	b = D(1,:)';
