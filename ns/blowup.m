function ifblow = blowup(u,v,p,s)

um = norm(u,inf);
vm = norm(v,inf);
pm = norm(p,inf);
sm = norm(s,inf);

ifblow = 0;

if    (um > 1e15); ifblow = 1;
elseif(vm > 1e15); ifblow = 1;
elseif(pm > 1e15); ifblow = 1;
elseif(sm > 1e15); ifblow = 1; end;

if(ifblow)
	['Blowup']
	[um,vm,pm,sm]
end

end
