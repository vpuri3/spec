%
ifvel  = 1;    % evolve velocity field per Navier-Stokes
ifpres = 1;    % project velocity field onto a div-free subspace
ifps   = 0;    % evolve passive scalar per advection diffusion eqn

Ex  = 4;
Ey  = 4;
nx1 = 8;
ny1 = nx1;
nx2 = nx1-2;
ny2 = ny1-2;
nxd = ceil(1.5*nx1)+rem(0.5*nx1,2);
nyd = ceil(1.5*ny1)+rem(0.5*ny1,2);

