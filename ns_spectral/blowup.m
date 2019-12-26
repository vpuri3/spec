function ifblow = blowup(u,v,p,s,Bv,Bp)

um = L2(u,Bv);
vm = L2(v,Bv);
pm = L2(p,Bp);
sm = L2(s,Bv);

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
