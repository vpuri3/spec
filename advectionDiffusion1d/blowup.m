function ifblow = blowup(u)

um = dot(u,u);

ifblow = 0;

if (um > 1e20); ifblow = 1; end;

if(ifblow)
	['Blowup'],[um]
end

end
