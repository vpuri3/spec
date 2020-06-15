function omega = vort(ux,uy,Ir,Is,Dr,Ds,rx,ry,sx,sy)

	[uxdx,uxdy] = grad(ux,Ir,Is,Dr,Ds,rx,ry,sx,sy);
	[uydx,uydy] = grad(uy,Ir,Is,Dr,Ds,rx,ry,sx,sy);

	omega = uxdx .* uydy - uxdy .* uydx;
