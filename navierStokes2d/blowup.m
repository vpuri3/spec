function ifblow = blowup(u,v,p,s)

um = dot(u,u);
vm = dot(v,v);
pm = dot(p,p);
sm = dot(s,s);

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
