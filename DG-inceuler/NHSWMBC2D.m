function Q = NHSWMBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
  
%  function [Q] = ForwardStepBC2D(xin, yin, nxin, nyin, mapI, mapO, mapW, mapC, Q, time);
% Purpose: Impose channel boundary conditions on 2D Euler equations on weak form


% extract conserved variables
h = Q(:,:,1); hu = Q(:,:,2); hv = Q(:,:,3); 

% Inflow conditions -- uniform inflow
%rhoin = gamma; uin = 3.0; vin = 0.0; pin = 1.0;
%Ein = pin/(gamma-1.0) + 0.5*rhoin*(uin^2+vin^2);

%rho(mapI) = rhoin; rhou(mapI) = rhoin*uin; rhov(mapI) = rhoin*vin; Ener(mapI) = Ein;

% Outflow conditions -- supersonic outflow ( do nothing )

% Wall conditions -- reflective, isothermal, i.e., n.u=0, T=T(t=0)
%hW = h(mapW); 
huW = hu(mapW); hvW = hv(mapW); 
nxW = nxin(mapW);   nyW = nyin(mapW);

% reverse flow in normal direction in ghost elements

%think you have to let the pressure step enforce no normal flow. for
%incompressible euler
%hu(mapW) = huW - 2*nxW.*(nxW.*huW + nyW.*hvW);
%hv(mapW) = hvW - 2*nyW.*(nxW.*huW + nyW.*hvW);
%hu(mapW) = huW;

%Seems to be missing the rho+ = rho- condition on the density.
%and the E+ = E- condition on the energy. See if you can reconcile this. --Derek.
%Reconciling: those types of conditions (zero neumann), require no work
%since at boundary nodes, vmapM = vmapP by default.


% pack modified conserved variables
Q(:,:,1) = h; Q(:,:,2) = hu; Q(:,:,3) = hv; %Q(:,:,4) = Ener;
return
