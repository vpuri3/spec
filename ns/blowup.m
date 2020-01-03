function ifblow = blowup(u,v,p,s,Bv,Bp,Qx1,Qy1,Qx2,Qy2)

um = L2(u,Bv,Qx1,Qy1);
vm = L2(v,Bv,Qx1,Qy1);
pm = L2(p,Bp,Qx2,Qy2);
sm = L2(s,Bv,Qx1,Qy1);

ifblow = 0;

if    (um > 1e20); ifblow = 1;
elseif(vm > 1e20); ifblow = 1;
elseif(pm > 1e20); ifblow = 1;
elseif(sm > 1e20); ifblow = 1; end;

if(ifblow)
	['Blowup']
	[um,vm,pm,sm]
end

end
