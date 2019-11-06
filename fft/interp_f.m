
function [um] = interp_f(m,un);

n  = size(un,1);

uh = rfft(un);
uh = [ uh; zeros(m-n,1) ];

um = irfft(uh);


