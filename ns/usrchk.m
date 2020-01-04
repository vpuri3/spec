%
if(blowup(vx,vy,pr,ps,Bm1,Bm2,Qx1,Qy1,Qx2,Qy2));it, return; end;

if(mod(it,5e1)==0 | time>=T-1e-6)
	% log
	%[it,L2(vx,Bm1),L2(vy,Bm1),L2(pr,Bm2),L2(ps,Bm1)]

	['infty kovazny normalized v-ve']
	[max(max(abs(vx-vxe))),max(max(abs(vy-vye)))] ./...
	[max(max(abs(vxe))),max(max(abs(vye)))]

	% vis
	om = vort(vx,vy,Qx1,Qy1,Dxm1,Dym1,rxm1,rym1,sxm1,sym1);

	contour(xm1g,ym1g,vx,20); view(2);colorbar
	%mesh(xm1g,ym1g,ps); view(3); if(it==nt) max(max(abs(ps-psb)));end
   	title([casename,', t=',num2str(time,'%4.2f'),' i=',num2str(it)]);
	drawnow
	if(T~=0) mov = [mov,getframe(fig)]; end
end

if(it == nt)
	['Finished Timestepping']
	%['Energy in vx,vy,pr,ps'],[L2(vx,Bm1),L2(vy,Bm1),L2(pr,Bm2),L2(ps,Bm1)]
	
	% play movie
	%movie(fig,mov,-2,40);
	
	% save gif
	
	%gname = [cname,'.gif'];
	%fps   = 40;
	%mov   = [mov,flip(mov)];
	%
	%for i=1:length(mov)
	%	f = mov(i);
	%	[img,cmap] = rgb2ind(f.cdata,256);
	%	if i==1 imwrite(img,cmap,gname,'gif','DelayTime',1/fps,'LoopCount',Inf)
	%	else imwrite(img,cmap,gname,'gif','WriteMode','append','DelayTime',1/fps)
	%	end
	%end
end
