function ifblow = blowup(u,v,p)

um = dot(u,u);
vm = dot(v,v);
pm = dot(p,p);

ifblow = 0;

if(um > 1e20); ifblow = 1; end;
if(vm > 1e20); ifblow = 1; end;
if(pm > 1e20); ifblow = 1; end;

if(ifblow)
	['Blowup']
	[um,vm,pm]
end

end
