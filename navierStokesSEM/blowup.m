function blowup(u,v,p,s,istep)

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
	format compact
	msg=['Step ',num2str(istep),', u v p ps',newline,num2str([um,vm,pm,sm])];
	[um,vm,pm,sm];
	error(msg);
end

end
