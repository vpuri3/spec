function omega = vort(ux,uy,Qx,Qy,Dr,Ds,rx,ry,sx,sy)

	uux = ABu(Qy,Qx,ux);
	uuy = ABu(Qy,Qx,uy);

	[uxdx,uxdy] = grad(uux,Dr,Ds,rx,ry,sx,sy);
	[uydx,uydy] = grad(uuy,Dr,Ds,rx,ry,sx,sy);

	omega = uxdx .* uydy - uxdy .* uydx;

	omega = ABu(Qy',Qx',omega);
