function c_grad_uh = advec_op(uh,Cx,Cy,Bd,Dx,Dy,Jx,Jy)
[Mx,My] = size(Bd);
[My,Ny] = size(Jy); Nx=size(uh,1); pad=Mx-Nx; 
% Cx,Cy formed on dealising mesh
% computing derv. and interpolating to M pts
Dyuh =    uh*(Jy*Dy)'; if pad>0; Dyuh=[Dyuh; zeros(pad,My)]; end;
Dxuh = Dx*uh*Jy';      if pad>0; Dxuh=[Dxuh; zeros(pad,My)]; end;

c_grad_uh=Cy.*irfft(Dyuh)+Cx.*irfft(Dxuh);  %% Collocate in phys space
c_grad_uh=rfft(c_grad_uh);                  %% Phys-to-wave transform
c_grad_uh=Bd .* c_grad_uh;                  %% Collocation w mass matrix
c_grad_uh=c_grad_uh*Jy;                     %% Interpolating to N pts (y)
c_grad_uh=c_grad_uh(1:Nx,:);                %% Interpolating to N pts (x)
