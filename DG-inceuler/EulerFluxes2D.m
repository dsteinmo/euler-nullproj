function [F,G,rho,u,v] = EulerFluxes2D(Q)
     
% function [F,G,rho,u,v,p] = EulerFluxes2D(Q, gamma)
% Purpose: evaluate primitive variables and Euler flux functions

% extract conserved variables
rho = Q(:,:,1); u = Q(:,:,2); v = Q(:,:,3); 
%Ener = Q(:,:,4);

%p = (gamma-1)*(Ener - 0.5*(rhou.*u + rhov.*v)); %no thermodynamic pressure

% compute flux functions
F = zeros(size(Q)); 
F(:,:,1) = rho.*u;
F(:,:,2) = u.*u;
F(:,:,3) = u.*v;
F(:,:,4) = zeros(size(u));

G = zeros(size(Q));
G(:,:,1) = rho.*v;
G(:,:,2) = u.*v; 
G(:,:,3) = v.*v; 
G(:,:,4) = zeros(size(u));
return;
