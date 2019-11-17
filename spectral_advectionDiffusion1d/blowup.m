function ifblow = blowup(u,s)

um = dot(u,u);
sm = dot(s,s);

ifblow = 0;

if    (um > 1e20); ifblow = 1;
elseif(sm > 1e20); ifblow = 1; end;

if(ifblow)
	['Blowup']
	[um,sm]
end

end
