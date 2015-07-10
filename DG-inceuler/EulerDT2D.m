%function dt = EulerDT2D(Q, gamma)
function dt = EulerDT2D(Q,g,H)

% function dt = EulerDT2D(Q, gamma)
% purpose: compute the time step dt for the compressible Euler equations

Globals2D;

rho = Q(:,:,1); 
u = Q(:,:,2); v = Q(:,:,3); 


u = u(vmapM); v = v(vmapM); 

drho = max(abs(rho(:)))-min(abs(rho(:)));
Nscale = sqrt(g*drho/H);
c0_bound = Nscale*H;

dt = 1/max( ((N+1)^2)*.5*abs(Fscale(:)).*(sqrt ( c0_bound+ u(:).^2 + v(:).^2 ) ));

%or use continuous version... (a bit more involved)
% rhor = Dr*rho; 
% rhos = Ds*rho;
% rhoy = ry.*rhor + sy.*rhos;
% Nfreq = sqrt(-g*rhoy);
% 
% Nfreq = Nfreq(vmapM);
% 
% c0_bound = Nfreq*H;
% 
% dt = 1/max( ((N+1)^2)*.5*Fscale(:).*(sqrt ( c0_bound(:)+ u(:).^2 + v(:).^2 ) ))
return
