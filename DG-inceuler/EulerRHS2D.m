function [rhsQ] = EulerRHS2D(Q,time, ExactSolutionBC,g)

% function [rhsQ] = EulerRHS2D(Q,time, ExactSolutionBC);
% Purpose: Evaluate RHS in 2D Euler equations, discretized on weak form
%             with a local Lax-Friedrich flux

Globals2D;

vmapM = reshape(vmapM, Nfp*Nfaces, K); vmapP = reshape(vmapP, Nfp*Nfaces, K);

% 1. Compute volume contributions (NOW INDEPENDENT OF SURFACE TERMS)
[F,G,rho,u,v] = EulerFluxes2D(Q);

% Compute weak derivatives
rhsQ = zeros(size(Q));
for n=1:3
  dFdr = Drw*F(:,:,n); dFds = Dsw*F(:,:,n);
  dGdr = Drw*G(:,:,n); dGds = Dsw*G(:,:,n);
  rhsQ(:,:,n) = (rx.*dFdr + sx.*dFds) + (ry.*dGdr + sy.*dGds);
end
    
% 2. Compute surface contributions 
% 2.1 evaluate '-' and '+' traces of conservative variables
%QM = zeros(size(Q));
%QP = zeros(size(Q));
for n=1:3
  Qn = Q(:,:,n);
  QM(:,:,n) = Qn(vmapM); QP(:,:,n) = Qn(vmapP);
end

% 2.2 set boundary conditions by modifying positive traces
if(~isempty(ExactSolutionBC))
  QP = feval(ExactSolutionBC, Fx, Fy, nx, ny, mapI, mapO, mapW, mapC, QP, time);
end

% 2.3 evaluate primitive variables & flux functions at '-' and '+' traces
[fM,gM,rhoM,uM,vM] = EulerFluxes2D(QM);
[fP,gP,rhoP,uP,vP] = EulerFluxes2D(QP);

%c0=1e3;

% 2.4 Compute local Lax-Friedrichs/Rusonov numerical fluxes
lambda = max( sqrt(uM.^2+vM.^2) ,  ...
	      sqrt(uP.^2+vP.^2));
%lambda = max( sqrt(uM.^2+vM.^2) + sqrt(abs(gamma*pM./rhoM)),  ...
%	      sqrt(uP.^2+vP.^2) + sqrt(abs(gamma*pP./rhoP)));
lambda = reshape(lambda, Nfp, Nfaces*K);
lambda = ones(Nfp, 1)*max(lambda, [], 1); 
lambda = reshape(lambda, Nfp*Nfaces, K);

%lambda = 0*lambda;

% 2.5 Lift fluxes
%lambda =0; 
for n=1:3 
  nflux = nx.*(fP(:,:,n) + fM(:,:,n)) + ny.*(gP(:,:,n) + gM(:,:,n)) + ...
      lambda.*(QM(:,:,n) - QP(:,:,n));
  rhsQ(:,:,n) = rhsQ(:,:,n) - LIFT*(Fscale.*nflux/2);
end

%add Dunphy's source term (makes code not blow-up for active tracers).
% dudr = Drw*Q(:,:,2); duds = Dsw*Q(:,:,2);
% dvdr = Drw*Q(:,:,3); dvds = Dsw*Q(:,:,3);
% mydiv = (rx.*dudr + sx.*duds) + (ry.*dvdr + sy.*dvds);

%above line seems to work

%SRCT = Q(:,:,1).*(mydiv - LIFT*(Fscale.*(nx.*(uM+uP) + ny.*(vM+vP))/2));
%this line does not work
%This doesn't look like it has anything to do with skew-symmetric forms.
%I think del . u = 0 isn't being dealt with properly across interfaces,
%so adding this term in effectively un-does that.


% u=Q(:,:,2);v=Q(:,:,3);
% udotn = u(vmapM).*nx+v(vmapM).*ny;
% 
%cancel spurious terms of the form phi*(div u)
%think doing this violates entropy
% for n=1:3
%      rhsQ(:,:,n) = rhsQ(:,:,n) - Q(:,:,n).*(mydiv - LIFT*(Fscale.*(nx.*(uM+uP) + ny.*(vM+vP))/2));
%      %rhsQ(:,:,n) = rhsQ(:,:,n) - Q(:,:,n).*mydiv; %results in insta-death
%      %rhsQ(:,:,n) = rhsQ(:,:,n) - Q(:,:,n).*(mydiv - ...
%      %             LIFT*(Fscale.*(nx.*((udotn>=0).*uM + (udotn<0).*uP))) ...
%      %                           +ny.*((udotn>=0).*vM + (udotn<0).*vP)); %dies horribly
%      
% end

%add buoyancy source term - Not necessary if hydrostatic pressure 
%handled separately with a z-integration.
rhsQ(:,:,3) = rhsQ(:,:,3) - g*rho;

%'hybridize' bit
%rotate to tangential/normal
%note: can't rotate whole RHS like this, can only do it with the traces
%algo: get traces, rotate, dss, rotate traces back, replace appropriate parts of
%whole RHS by traces
rhsrho =rhsQ(:,:,1);
rhsu= rhsQ(:,:,2);
rhsv= rhsQ(:,:,3);



method='dg';
%method='hybrid';
%method='specel';

%full spectral elements
if strcmp(method,'specel')
    rhsrho(vmapM) = .5*(rhsrho(vmapM)+rhsrho(vmapP));
    rhsrho(vmapP) = rhsrho(vmapM);

    rhsu(vmapM) = .5*(rhsu(vmapM)+rhsu(vmapP));
    rhsu(vmapP) = rhsu(vmapM);

    rhsv(vmapM) = .5*(rhsv(vmapM)+rhsv(vmapP));
    rhsv(vmapP) = rhsv(vmapM);
elseif strcmp(method,'hybrid')
    rhsuM = rhsu(vmapM); rhsvM = rhsv(vmapM);
    rhsuP = rhsu(vmapP); rhsvP = rhsv(vmapP);
    rhsunM = nx.*rhsuM + ny.*rhsvM; 
    rhsutM =-ny.*rhsuM + nx.*rhsvM; 
    rhsunP = nx.*rhsuP + ny.*rhsvP; 
    rhsutP =-ny.*rhsuP + nx.*rhsvP; 
    
    %avg. in normal direction (replace with upwind?)
    %rhsunM = .5*(rhsunM + rhsunP);
    %rhsunP = rhsunM;
    
    %upwind
    u=Q(:,:,2);v=Q(:,:,3);
    udotn = u(vmapM).*nx+v(vmapM).*ny;
    rhsunM = (udotn>=0).*rhsunM + (udotn<0).*rhsunP;
    rhsunP = rhsunM;
    
    
    rhsuM = nx.*rhsunM - ny.*rhsutM;
    rhsvM = ny.*rhsunM + nx.*rhsutM;
    rhsuP = nx.*rhsunP - ny.*rhsutP;
    rhsvP = ny.*rhsunP + nx.*rhsutP;
    
    rhsu(vmapM) = rhsuM;
    rhsu(vmapP) = rhsuP;
    rhsv(vmapM) = rhsvM;
    rhsv(vmapP) = rhsvP;
    
end

rhsQ(:,:,1) = rhsrho;
rhsQ(:,:,2) = rhsu;
rhsQ(:,:,3) = rhsv;


return;


