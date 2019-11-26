%
function [Hu] =  hmhltz(u,visc,b0,Bd,Jr,Js,Dr,Ds,G11,G12,G22)
	Hu =      visc*lapl(u,Jr,Js,Dr,Ds,G11,G12,G22);
	Hu = Hu +   b0*mass(u,Bd,Jr,Js);
end
