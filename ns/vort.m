function omega = vort(ux,uy,Dr,Ds,rx,ry,sx,sy)

	[uxdx,uxdy] = grad(ux,Dr,Ds,rx,ry,sx,sy);
	[uydx,uydy] = grad(uy,Dr,Ds,rx,ry,sx,sy);

	omega = uxdx .* uydy - uxdy .* uydx;
